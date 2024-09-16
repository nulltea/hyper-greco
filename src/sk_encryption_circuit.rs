use crate::constants::BfvSkEncryptConstans;
use crate::lasso::LassoPreprocessing;
use crate::{
    lasso::{table::range::RangeLookup, LassoNode},
    poly::Poly,
    transcript::Keccak256Transcript,
};
use gkr::izip_eq;
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
use itertools::chain;
use plonkish_backend::pcs::PolynomialCommitmentScheme;
use plonkish_backend::poly::multilinear::MultilinearPolynomial;
use plonkish_backend::util::hash::{Keccak256, Output};
use rand::RngCore;
use rayon::iter::ParallelIterator;
use serde::Deserialize;
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

/// `BfvSkEncryptionCircuit` is a circuit that checks the correct formation of a ciphertext resulting from BFV secret key encryption
/// All the polynomials coefficients and scalars are normalized to be in the range `[0, p)` where p is the modulus of the prime field of the circuit
///
/// # Parameters:
/// * `s`: secret polynomial, sampled from ternary distribution.
/// * `e`: error polynomial, sampled from discrete Gaussian distribution.
/// * `k1`: scaled message polynomial.
/// * `r2is`: list of r2i polynomials for each i-th CRT basis .
/// * `r1is`: list of r1i polynomials for each CRT i-th CRT basis.
/// * `ais`: list of ai polynomials for each CRT i-th CRT basis.
/// * `ct0is`: list of ct0i (first component of the ciphertext cti) polynomials for each CRT i-th CRT basis.
#[derive(Deserialize, Clone)]
pub struct BfvSkEncryptArgs {
    s: Vec<String>,
    e: Vec<String>,
    k1: Vec<String>,
    r2is: Vec<Vec<String>>,
    r1is: Vec<Vec<String>>,
    ais: Vec<Vec<String>>,
    ct0is: Vec<Vec<String>>,
}

pub struct BfvEncryptBlock<Params: BfvSkEncryptConstans<K>, const K: usize> {
    num_reps: usize,
    _marker: std::marker::PhantomData<Params>,
}

impl<Params: BfvSkEncryptConstans<K>, const K: usize> BfvEncryptBlock<Params, K> {
    pub const fn log2_size(&self) -> usize {
        Params::N_LOG2 + 1
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
        let poly_log2_size = Params::N_LOG2;
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
                        relay_mul_const((0, j), F::from_str_vartime(Params::K0IS[i]).unwrap())
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
                        relay_mul_const((i, j), F::from_str_vartime(Params::QIS[i]).unwrap())
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
                Params::R1_BOUNDS.iter().take(self.num_reps).copied(),
                iter::repeat(Params::R2_BOUNDS[0]).take(r2is_chunks.len()),
                vec![Params::S_BOUND, Params::E_BOUND, Params::K1_BOUND],
            ]
            .enumerate()
            .flat_map(|(i, bound)| {
                (0..(1usize << log2_size)).map(move |j| relay_add_const((i, j), F::from(bound)))
            })
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
                Params::R1_BOUNDS
                    .iter()
                    .take(self.num_reps)
                    .flat_map(|&bound| iter::repeat(RangeLookup::id_for(bound * 2 + 1))
                        .take(1 << log2_size)),
                Params::R2_BOUNDS
                    .iter()
                    .take(self.num_reps)
                    .flat_map(|&bound| iter::repeat(RangeLookup::id_for(bound * 2 + 1))
                        .take(1 << r2i_log2_size)),
                iter::repeat(RangeLookup::id_for(Params::S_BOUND * 2 + 1)).take(1 << log2_size),
                iter::repeat(RangeLookup::id_for(Params::E_BOUND * 2 + 1)).take(1 << log2_size),
                iter::repeat(RangeLookup::id_for(Params::K1_BOUND * 2 + 1)).take(1 << log2_size),
            ]
            .collect_vec();
            let num_vars = lookups.len().next_power_of_two().ilog2() as usize;
            circuit.insert(LassoNode::<F, E, C, M>::new(
                preprocessing,
                num_vars,
                lookups,
            ))
        };
        r1is.iter()
            .take(self.num_reps)
            .for_each(|&r1i| circuit.connect(r1i, lasso_inputs_batched));

        r2is_chunks
            .iter()
            .for_each(|&r2i| circuit.connect(r2i, lasso_inputs_batched));

        connect!(circuit {
            lasso_inputs_batched <- s, e, k1;
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

pub struct BfvEncrypt<Params: BfvSkEncryptConstans<K>, const K: usize> {
    block: BfvEncryptBlock<Params, K>,
}

impl<Params: BfvSkEncryptConstans<K>, const K: usize> BfvEncrypt<Params, K> {
    pub fn new(num_reps: usize) -> Self {
        Self {
            block: BfvEncryptBlock {
                num_reps,
                _marker: std::marker::PhantomData,
            },
        }
    }

    pub const fn log2_size(&self) -> usize {
        Params::N_LOG2 + 1
    }

    #[allow(clippy::type_complexity)]
    pub fn setup<
        F: PrimeField,
        E: ExtensionField<F>,
        Pcs: PolynomialCommitmentScheme<F, Polynomial = MultilinearPolynomial<F>>,
    >(
        &self,
        _rng: impl RngCore + Clone,
    ) -> (ProverKey<F, E>, VerifierKey<F, E>) {
        let mut lasso_preprocessing = LassoPreprocessing::<F, E>::preprocess::<C, M>(chain![
            [
                RangeLookup::new_boxed(Params::S_BOUND * 2 + 1),
                RangeLookup::new_boxed(Params::E_BOUND * 2 + 1),
                RangeLookup::new_boxed(Params::K1_BOUND * 2 + 1)
            ],
            Params::R1_BOUNDS
                .iter()
                .take(self.block.num_reps)
                .map(|&bound| RangeLookup::new_boxed(bound * 2 + 1)),
            Params::R2_BOUNDS
                .iter()
                .take(self.block.num_reps)
                .map(|&bound| RangeLookup::new_boxed(bound * 2 + 1))
        ]);

        let lasso_verifier = lasso_preprocessing.to_verifier_preprocessing();

        let pk = lasso_preprocessing;
        let vk = lasso_verifier;

        (pk, vk)
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

    pub fn prove<
        F: PrimeField,
        Pcs: PolynomialCommitmentScheme<
            F,
            Polynomial = MultilinearPolynomial<F>,
            CommitmentChunk = Output<Keccak256>,
        >,
    >(
        &self,
        args: &BfvSkEncryptArgs,
        pk: ProverKey<F, F>,
    ) -> Vec<u8> {
        let preprocessing = pk;
        let mut transcript = Keccak256Transcript::<Vec<u8>>::default();

        let circuit = info_span!("init circuit").in_scope(|| {
            let mut circuit = Circuit::<F, F>::default();
            self.configure(&mut circuit, preprocessing);
            circuit
        });

        let (values, output_claims) = info_span!("wintess gen").in_scope(|| {
            let (inputs, ctis_poly) =
                info_span!("parse inputs").in_scope(|| self.get_inputs::<F, F>(args));

            let values = info_span!("eval circuit").in_scope(|| circuit.evaluate(inputs));

            let ct0is_claim = info_span!("eval output").in_scope(|| {
                let point = transcript.squeeze_challenges(self.ct0is_log2_size());
                let value = ctis_poly.evaluate(&point);
                EvalClaim::new(point.clone(), value)
            });

            let output_claims = vec![EvalClaim::new(vec![], F::ZERO), ct0is_claim];

            (values, output_claims)
        });

        let _claims = info_span!("GKR prove")
            .in_scope(|| gkr::prove_gkr(&circuit, &values, &output_claims, &mut transcript))
            .unwrap();

        transcript.into_proof()
    }

    pub fn verify<
        F: PrimeField,
        Pcs: PolynomialCommitmentScheme<
            F,
            Polynomial = MultilinearPolynomial<F>,
            CommitmentChunk = Output<Keccak256>,
        >,
    >(
        &self,
        vk: VerifierKey<F, F>,
        inputs: Vec<BoxMultilinearPoly<'static, F, F>>,
        proof: &[u8],
        ct0is: Vec<Vec<String>>,
    ) {
        let preprocessing = vk;
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

            vec![EvalClaim::new(vec![], F::ZERO), ct0is_claim]
        });

        let circuit = info_span!("init circuit").in_scope(|| {
            let mut circuit = Circuit::<F, F>::default();
            self.configure(&mut circuit, preprocessing);
            circuit
        });

        let input_claims = info_span!("GKR verify")
            .in_scope(|| verify_gkr(&circuit, &output_claims, &mut transcript).unwrap());

        izip_eq!(inputs, input_claims).for_each(|(input, claims)| {
            claims
                .iter()
                .for_each(|claim| assert_eq!(input.evaluate(claim.point()), claim.value()))
        });
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
    use crate::generate_sk_enc_test;

    use super::*;
    use gkr::util::dev::seeded_std_rng;
    use goldilocks::Goldilocks;

    use paste::paste;
    use plonkish_backend::{pcs::multilinear::MultilinearBrakedown, util::code::BrakedownSpec6};
    use std::{fs::File, io::Read};
    use tracing::info_span;
    use tracing_forest::ForestLayer;
    use tracing_subscriber::{layer::SubscriberExt, EnvFilter, Registry};

    pub type Brakedown<F> =
        MultilinearBrakedown<F, plonkish_backend::util::hash::Keccak256, BrakedownSpec6>;

    generate_sk_enc_test!(goldilocks, Goldilocks, Brakedown<Goldilocks>, 1024, 1, 27);

    generate_sk_enc_test!(goldilocks, Goldilocks, Brakedown<Goldilocks>, 2048, 1, 52);

    generate_sk_enc_test!(goldilocks, Goldilocks, Brakedown<Goldilocks>, 4096, 2, 55);

    generate_sk_enc_test!(goldilocks, Goldilocks, Brakedown<Goldilocks>, 8192, 4, 55);

    generate_sk_enc_test!(goldilocks, Goldilocks, Brakedown<Goldilocks>, 16384, 8, 54);
}
