use crate::{poly::Poly, transcript::Keccak256Transcript};
use bfv::BfvParameters;
use gkr::izip_eq;
use gkr::util::arithmetic::radix2_fft;
use gkr::{
    chain_par,
    circuit::{
        connect,
        node::{EvalClaim, FftNode, InputNode, VanillaGate, VanillaNode},
        Circuit, NodeId,
    },
    ff_ext::ff::PrimeField,
    poly::{box_dense_poly, BoxMultilinearPoly},
    transcript::Transcript,
    util::{arithmetic::ExtensionField, Itertools},
    verify_gkr,
};
use itertools::{chain, izip};
use lasso_gkr::{table::range::RangeLookup, LassoNode, LassoPreprocessing};
use num_bigint::BigInt;
use num_traits::ToPrimitive;
use plonkish_backend::util::arithmetic::root_of_unity;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use std::cmp::min;
use std::iter;
use tracing::info_span;

const LIMB_BITS: usize = 16;
const C: usize = 4;
const M: usize = 1 << LIMB_BITS;

pub type ProverKey<
    F,
    E,
    // Pcs: PolynomialCommitmentScheme<
    //     F,
    //     Polynomial = MultilinearPolynomial<F>,
    //     CommitmentChunk = Output<Keccak256>,
    // >,
> = LassoPreprocessing<F, E>;

pub type VerifierKey<
    F,
    E,
    // Pcs: PolynomialCommitmentScheme<
    //     F,
    //     Polynomial = MultilinearPolynomial<F>,
    //     CommitmentChunk = Output<Keccak256>,
    // >,
> = LassoPreprocessing<F, E>;

pub type BfvSkEncryptArgs = bfv_witgen::InputValidationVectors;

pub struct BfvEncryptBlock {
    params: BfvParameters,
    bounds: bfv_witgen::InputValidationBounds,
    num_reps: usize,
}

impl BfvEncryptBlock {
    pub const fn log2_size(&self) -> usize {
        self.params.degree.ilog2() as usize + 1
    }

    // single block
    pub fn configure<F: PrimeField, E: ExtensionField<F>>(
        &self,
        circuit: &mut Circuit<F, E>,
        s: NodeId,
        e: NodeId,
        k1: NodeId,
        preprocessing: LassoPreprocessing<F, E>,
    ) -> NodeId {
        let poly_log2_size = self.params.degree.ilog2() as usize;
        let log2_size = self.log2_size();

        let es = {
            let gates = (0..self.num_reps)
                .flat_map(|_| (0..(1usize << log2_size)).map(move |j| VanillaGate::relay((0, j))))
                .collect_vec();

            circuit.insert(VanillaNode::new(1, log2_size, gates.clone(), 1))
        };

        let k1kis = {
            let gates = (0..self.num_reps)
                .flat_map(|i| {
                    (0..(1usize << log2_size)).map(move |j| {
                        relay_mul_const(
                            (0, j),
                            F::from_u128(self.bounds.k01s[i].to_u128().unwrap()),
                        )
                    })
                })
                .collect_vec();

            circuit.insert(VanillaNode::new(1, log2_size, gates.clone(), 1))
        };

        connect!(circuit {
            es <- e;
            k1kis <- k1;
        });

        let ais = iter::repeat_with(|| circuit.insert(InputNode::new(log2_size, 1)))
            .take(self.num_reps)
            .collect_vec();

        let r1is = iter::repeat_with(|| circuit.insert(InputNode::new(log2_size, 1)))
            .take(self.num_reps)
            .collect_vec();

        let r1iqis = {
            let r1i_size = 1usize << log2_size;
            let gates = (0..self.num_reps)
                .flat_map(|i| {
                    (0..r1i_size).map(move |j| {
                        relay_mul_const((i, j), F::from_u128(self.bounds.qis[i].to_u128().unwrap()))
                    })
                })
                .collect_vec();

            circuit.insert(VanillaNode::new(self.num_reps, log2_size, gates.clone(), 1))
        };

        r1is.iter()
            .take(self.num_reps)
            .for_each(|&r1i| circuit.connect(r1i, r1iqis));

        let r2is = circuit.insert(InputNode::new(poly_log2_size, self.num_reps));

        let r2is_log2_sise = self.log2_size_with_num_reps(poly_log2_size);
        let r2is_chunks = (0..1 << r2is_log2_sise)
            .chunks(1 << log2_size)
            .into_iter()
            .map(|chunk| {
                let mut gates = chunk.map(move |j| VanillaGate::relay((0, j))).collect_vec();
                gates.resize(1 << log2_size, VanillaGate::constant(F::ZERO));

                let node = circuit.insert(VanillaNode::new(1, r2is_log2_sise, gates, 1));
                circuit.connect(r2is, node);
                node
            })
            .collect_vec();

        let lasso_inputs_batched = {
            let gates = chain![
                vec![self.bounds.s, self.bounds.e, self.bounds.k1],
                self.bounds.r1.iter().take(self.num_reps).copied(),
            ]
            .enumerate()
            .flat_map(|(i, bound)| {
                (0..(1usize << log2_size)).map(move |j| relay_add_const((i, j), F::from(bound)))
            })
            .chain(
                self.bounds
                    .r2
                    .iter()
                    .take(self.num_reps)
                    .copied()
                    .enumerate()
                    .flat_map(|(i, bound)| {
                        (0..(1usize << poly_log2_size)).map(move |j| {
                            let input_gate = self.num_reps + 3 + i / 2; // (s,e,k1) + r1s + r2s chunks (each 2x of r2i)
                            relay_add_const((input_gate, j), F::from(bound))
                        })
                    }),
            )
            .collect_vec();

            circuit.insert(VanillaNode::new(
                r2is_chunks.len() + self.num_reps + 3,
                log2_size,
                gates,
                1,
            ))
        };
        let lasso_ranges = {
            let r2i_log2_size = if self.num_reps == 1 {
                log2_size // since zero-padded to log2_size in r2is_chunks
            } else {
                poly_log2_size
            };
            let lookups = chain![
                iter::repeat(RangeLookup::id_for(self.bounds.s * 2 + 1)).take(1 << log2_size),
                iter::repeat(RangeLookup::id_for(self.bounds.e * 2 + 1)).take(1 << log2_size),
                iter::repeat(RangeLookup::id_for(self.bounds.k1 * 2 + 1)).take(1 << log2_size),
                self.bounds
                    .r1
                    .iter()
                    .take(self.num_reps)
                    .flat_map(|&bound| iter::repeat(RangeLookup::id_for(bound * 2 + 1))
                        .take(1 << log2_size)),
                self.bounds
                    .r2
                    .iter()
                    .take(self.num_reps)
                    .flat_map(|&bound| iter::repeat(RangeLookup::id_for(bound * 2 + 1))
                        .take(1 << r2i_log2_size)),
            ]
            .collect_vec();
            let num_vars = lookups.len().next_power_of_two().ilog2() as usize;
            circuit.insert(LassoNode::<F, E, C, M>::new(
                preprocessing,
                num_vars,
                lookups,
            ))
        };

        connect!(circuit {
            lasso_inputs_batched <- s, e, k1;
        });

        r1is.iter()
            .take(self.num_reps)
            .for_each(|&r1i| circuit.connect(r1i, lasso_inputs_batched));

        r2is_chunks
            .iter()
            .for_each(|&r2i| circuit.connect(r2i, lasso_inputs_batched));

        connect!(circuit {
            lasso_ranges <- lasso_inputs_batched;
        });

        let s_eval = circuit.insert(FftNode::forward(log2_size));
        circuit.connect(s, s_eval);

        let s_eval_copy = circuit.insert(VanillaNode::new(
            1,
            log2_size,
            (0..1usize << log2_size)
                .map(|i| VanillaGate::relay((0, i)))
                .collect_vec(),
            1,
        ));
        circuit.connect(s_eval, s_eval_copy);

        let sai_par = {
            let gates = (0..self.num_reps)
                .flat_map(|i| (0..(1usize << log2_size)).map(move |j| VanillaGate::relay((i, j))))
                .collect_vec();

            circuit.insert(VanillaNode::new(self.num_reps, log2_size, gates.clone(), 1))
        };

        for &ai in ais.iter().take(self.num_reps) {
            let gates = (0..1usize << log2_size)
                .map(|i| VanillaGate::mul((0, i), (1, i)))
                .collect_vec();
            let ai_eval = circuit.insert(FftNode::forward(log2_size));
            let sai_eval = circuit.insert(VanillaNode::new(2, log2_size, gates, 1));
            let sai = circuit.insert(FftNode::inverse(log2_size));

            connect!(circuit {
                ai_eval <- ai;
                sai_eval <- s_eval_copy, ai_eval;
                sai <- sai_eval;
            });

            circuit.connect(sai, sai_par);
        }

        let r2i_cyclo = {
            let r2i_size = (1usize << poly_log2_size) - 1;
            let gates = chain![
                (0..r2i_size).map(|i| VanillaGate::relay((0, i))),
                [VanillaGate::constant(F::ZERO)],
                (0..r2i_size).map(|i| VanillaGate::relay((0, i))),
                [VanillaGate::constant(F::ZERO)]
            ]
            .collect_vec();

            circuit.insert(VanillaNode::new(
                1,
                poly_log2_size,
                gates.clone(),
                self.num_reps,
            ))
        };

        let sum = {
            let gates = (0..1usize << log2_size)
                .map(|i| VanillaGate::sum(vec![(0, i), (1, i), (2, i), (3, i), (4, i)]))
                .collect();
            circuit.insert(VanillaNode::new(5, log2_size, gates, self.num_reps))
        };

        connect!(circuit {
            r2i_cyclo <- r2is;
            sum <- sai_par, es, k1kis, r1iqis, r2i_cyclo;
        });

        k1
    }

    fn log2_size_with_num_reps(&self, poly_log2_size: usize) -> usize {
        poly_log2_size + self.num_reps.ilog2() as usize
    }
}

pub struct BfvEncrypt {
    block: BfvEncryptBlock,
}

impl BfvEncrypt {
    pub fn new(
        params: BfvParameters,
        bounds: bfv_witgen::InputValidationBounds,
        num_reps: usize,
    ) -> Self {
        Self {
            block: BfvEncryptBlock {
                params,
                bounds,
                num_reps,
            },
        }
    }

    pub const fn log2_size(&self) -> usize {
        self.block.log2_size()
    }

    #[allow(clippy::type_complexity)]
    pub fn setup<
        F: PrimeField,
        E: ExtensionField<F>,
        // Pcs: PolynomialCommitmentScheme<F, Polynomial = MultilinearPolynomial<F>>,
    >(
        &self,
    ) -> ProverKey<F, E> {
        LassoPreprocessing::<F, E>::preprocess::<C, M>(chain![
            [
                RangeLookup::new_boxed(self.block.bounds.s * 2 + 1),
                RangeLookup::new_boxed(self.block.bounds.e * 2 + 1),
                RangeLookup::new_boxed(self.block.bounds.k1 * 2 + 1)
            ],
            self.block
                .bounds
                .r1
                .iter()
                .take(self.block.num_reps)
                .map(|&bound| RangeLookup::new_boxed(bound * 2 + 1)),
            self.block
                .bounds
                .r2
                .iter()
                .take(self.block.num_reps)
                .map(|&bound| RangeLookup::new_boxed(bound * 2 + 1))
        ])
    }

    pub fn configure<F: PrimeField, E: ExtensionField<F>>(
        &self,
        circuit: &mut Circuit<F, E>,
        preprocessing: LassoPreprocessing<F, E>,
    ) -> NodeId {
        let log2_size = self.log2_size();

        let s = circuit.insert(InputNode::new(log2_size, 1));
        let e = circuit.insert(InputNode::new(log2_size, 1));
        let k1 = circuit.insert(InputNode::new(log2_size, 1));

        self.block.configure(circuit, s, e, k1, preprocessing)
    }

    pub fn get_inputs<F: PrimeField, E: ExtensionField<F>>(
        &self,
        args: &BfvSkEncryptArgs,
    ) -> (
        Vec<BoxMultilinearPoly<'static, F, E>>,
        BoxMultilinearPoly<'static, F, E>,
    ) {
        let log2_size = self.log2_size();

        let s = Poly::<F>::new_padded(args.s.clone(), log2_size);
        let e = Poly::<F>::new_shifted(args.e.clone(), (1 << log2_size) - 1);
        let k1 = Poly::<F>::new_shifted(args.k1.clone(), (1 << log2_size) - 1);

        let mut r2is = vec![];
        let mut r1is = vec![];
        let mut ais = vec![];
        let mut ct0is = vec![];

        for z in 0..min(args.ct0is.len(), self.block.num_reps) {
            let r2i = Poly::<F>::new(args.r2is[z].clone());
            r2is.push(r2i.to_vec());

            let r1i = Poly::<F>::new_padded(args.r1is[z].clone(), log2_size);
            r1is.push(r1i.to_vec());

            let ai = Poly::<F>::new_padded(args.ais[z].clone(), log2_size);
            ais.push(ai.to_vec());

            let ct0i = Poly::<F>::new_shifted(args.ct0is[z].clone(), 1 << log2_size);
            let mut ct0i = ct0i.as_ref()[1..].to_vec();
            ct0i.push(F::ZERO);
            ct0is.extend(ct0i);
        }

        let r2is = r2is
            .into_iter()
            .take(self.block.num_reps)
            .flat_map(|mut r2i| {
                r2i.push(F::ZERO);
                r2i
            })
            .collect_vec();

        let inputs = chain_par![[s.to_vec(), e.to_vec(), k1.to_vec()], ais, r1is, [r2is],]
            .map(box_dense_poly)
            .collect();

        let output = box_dense_poly(ct0is);

        (inputs, output)
    }

    pub fn build_circuit<F: PrimeField, E: ExtensionField<F>>(
        &self,
        pk: ProverKey<F, E>,
    ) -> Circuit<F, E> {
        let mut circuit = Circuit::<F, E>::default();
        self.configure(&mut circuit, pk);
        circuit
    }

    pub fn prove<F: PrimeField, E: ExtensionField<F>>(
        &self,
        args: &BfvSkEncryptArgs,
        circuit: &Circuit<F, E>,
    ) -> Vec<u8> {
        let mut transcript = Keccak256Transcript::<Vec<u8>>::default();

        let (values, output_claims) = info_span!("wintess gen").in_scope(|| {
            let (inputs, ctis_poly) = info_span!("parse inputs").in_scope(|| self.get_inputs(args));

            let values = info_span!("eval circuit").in_scope(|| circuit.evaluate(inputs));

            let ct0is_claim = info_span!("eval output").in_scope(|| {
                let point = transcript.squeeze_challenges(self.ct0is_log2_size());
                let value = ctis_poly.evaluate(&point);
                EvalClaim::new(point.clone(), value)
            });

            let output_claims = vec![EvalClaim::new(vec![], E::ZERO), ct0is_claim];

            (values, output_claims)
        });

        let _claims = info_span!("GKR prove")
            .in_scope(|| gkr::prove_gkr(circuit, &values, &output_claims, &mut transcript))
            .unwrap();

        transcript.into_proof()
    }

    pub fn verify<F: PrimeField, E: ExtensionField<F>>(
        &self,
        // inputs: Vec<BoxMultilinearPoly<'static, F, E>>,
        circuit: &Circuit<F, E>,
        ct0is: Vec<Vec<BigInt>>,
        proof: &[u8],
    ) {
        let mut transcript = Keccak256Transcript::from_proof(proof);

        let output_claims = info_span!("eval output claim").in_scope(|| {
            let ct0is_claim = {
                let point = transcript.squeeze_challenges(self.ct0is_log2_size());
                let ct0is = box_dense_poly(
                    ct0is
                        .into_iter()
                        .take(self.block.num_reps)
                        .flat_map(|ct0i| {
                            let ct0i = Poly::<F>::new_shifted(ct0i, 1 << self.log2_size());
                            let mut ct0i = ct0i.as_ref()[1..].to_vec();
                            ct0i.push(F::ZERO);
                            ct0i
                        })
                        .collect_vec(),
                );
                let value = ct0is.evaluate(&point);

                EvalClaim::new(point, value)
            };

            vec![EvalClaim::new(vec![], E::ZERO), ct0is_claim]
        });

        let _input_claims = info_span!("GKR verify")
            .in_scope(|| verify_gkr(&circuit, &output_claims, &mut transcript).unwrap());

        // izip_eq!(inputs, input_claims).for_each(|(input, claims)| {
        //     claims
        //         .iter()
        //         .for_each(|claim| assert_eq!(input.evaluate(claim.point()), claim.value()))
        // });
    }

    fn ct0is_log2_size(&self) -> usize {
        assert!(self.block.num_reps.is_power_of_two());
        self.log2_size() + self.block.num_reps.next_power_of_two().ilog2() as usize
    }
}

fn relay_mul_const<F>(w: (usize, usize), c: F) -> VanillaGate<F> {
    VanillaGate::new(None, vec![(Some(c), w)], Vec::new())
}

fn relay_add_const<F>(w: (usize, usize), c: F) -> VanillaGate<F> {
    VanillaGate::new(Some(c), vec![(None, w)], Vec::new())
}

#[cfg(test)]
mod test {
    use super::*;
    use bfv::{BfvParameters, Encoding, Plaintext, SecretKey};
    use goldilocks::{Goldilocks, GoldilocksExt2};
    use num_traits::FromPrimitive;
    use rand::{distributions::Uniform, rngs::StdRng, thread_rng, Rng, SeedableRng};
    use tracing::info_span;
    use tracing_forest::ForestLayer;
    use tracing_subscriber::{layer::SubscriberExt, EnvFilter, Registry};
    use traits::TryEncodingWithParameters;

    fn random_bit_vector(hamming_weight: usize, size: usize) -> Vec<u32> {
        let mut rng = thread_rng();

        let mut bit_vector = vec![0; size];
        (0..hamming_weight).for_each(|_| {
            let sample_index = rng.sample(Uniform::new(0, size));
            bit_vector[sample_index] = 1;
        });

        bit_vector
    }

    #[test]
    fn test_sk_enc_valid() {
        let env_filter = EnvFilter::builder()
            .with_default_directive(tracing::Level::INFO.into())
            .from_env_lossy();

        let subscriber = Registry::default()
            .with(env_filter)
            .with(ForestLayer::default());

        let _ = tracing::subscriber::set_global_default(subscriber);

        let mut rng = thread_rng();

        const RING_SIZE: usize = 32;

        let mut params = BfvParameters::new_with_primes(
            vec![1032193, 1073692673],
            vec![995329, 1073668097],
            40961,
            RING_SIZE,
        );
        params.enable_hybrid_key_switching_with_prime(vec![61441]);
        params.enable_pke();

        let N: u64 = params.degree as u64;
        let bounds = bfv_witgen::witness_bounds(&params).unwrap();
        let bfv = BfvEncrypt::new(params.clone(), bounds, 1);
        let mut pk = info_span!("setup").in_scope(|| bfv.setup::<Goldilocks, GoldilocksExt2>());

        let circuit = info_span!("init circuit").in_scope(|| bfv.build_circuit(pk));

        let hamming_weight = 1000;
        let vector_size = RING_SIZE;

        let args = {
            let sk = SecretKey::random_with_params(&params, &mut rng);

            let m: Vec<_> = (0..N).collect_vec();
            let pt = Plaintext::encode(&m, &params, Encoding::default());

            let p = BigInt::from_u64(18446744069414584321).unwrap();

            bfv_witgen::encrypt_with_witness(&params, &sk, pt, &mut rng, &p)
                .unwrap()
                .1
        };

        let proof = info_span!("FHE_enc prove")
            .in_scope(|| bfv.prove::<Goldilocks, GoldilocksExt2>(&args, &circuit));

        info_span!("FHE_enc verify")
            .in_scope(|| bfv.verify::<Goldilocks, GoldilocksExt2>(&circuit, args.ct0is, &proof));
    }

    //     // Goldilocks prime tests

    //     generate_sk_enc_test!(
    //         "goldilocks",
    //         Goldilocks,
    //         GoldilocksExt2,
    //         Brakedown<Goldilocks>,
    //         1024,
    //         1,
    //         27
    //     );

    //     generate_sk_enc_test!(
    //         "goldilocks",
    //         Goldilocks,
    //         GoldilocksExt2,
    //         Brakedown<Goldilocks>,
    //         2048,
    //         1,
    //         52
    //     );

    //     generate_sk_enc_test!(
    //         "goldilocks",
    //         Goldilocks,
    //         GoldilocksExt2,
    //         Brakedown<Goldilocks>,
    //         4096,
    //         2,
    //         55
    //     );

    //     generate_sk_enc_test!(
    //         "goldilocks",
    //         Goldilocks,
    //         GoldilocksExt2,
    //         Brakedown<Goldilocks>,
    //         8192,
    //         4,
    //         55
    //     );

    //     generate_sk_enc_test!(
    //         "goldilocks",
    //         Goldilocks,
    //         GoldilocksExt2,
    //         Brakedown<Goldilocks>,
    //         16384,
    //         8,
    //         54
    //     );

    //     generate_sk_enc_test!(
    //         "goldilocks",
    //         Goldilocks,
    //         GoldilocksExt2,
    //         Brakedown<Goldilocks>,
    //         32768,
    //         16,
    //         59
    //     );

    //     // Bn254 prime tests

    //     generate_sk_enc_test!("bn254", Fr, Fr, Brakedown<Fr>, 1024, 1, 27);

    //     generate_sk_enc_test!("bn254", Fr, Fr, Brakedown<Fr>, 2048, 1, 52);

    //     generate_sk_enc_test!("bn254", Fr, Fr, Brakedown<Fr>, 4096, 2, 55);

    //     generate_sk_enc_test!("bn254", Fr, Fr, Brakedown<Fr>, 8192, 4, 55);

    //     generate_sk_enc_test!("bn254", Fr, Fr, Brakedown<Fr>, 16384, 8, 54);

    //     generate_sk_enc_test!("bn254", Fr, Fr, Brakedown<Fr>, 32768, 16, 59);
}
