use crate::constants::sk_enc_constants_1024_2x55_65537::{
    E_BOUND, K0IS, K1_BOUND, N, QIS, R1_BOUNDS, R2_BOUNDS,
};
use crate::pcs::{self, Ligero};
use crate::poly::Poly;
use gkr::izip_par;
use gkr::transcript::TranscriptWrite;
use gkr::util::arithmetic::radix2_fft;
use gkr::{
    chain_par,
    circuit::{
        connect,
        node::{FftNode, InputNode, LogUpNode, VanillaGate, VanillaNode},
        Circuit, NodeId,
    },
    ff_ext::ff::PrimeField,
    poly::{box_dense_poly, BoxMultilinearPoly},
    util::{
        arithmetic::{squares, ExtensionField},
        izip, Itertools,
    },
};
use itertools::chain;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use serde::Deserialize;
use std::cmp::{max, min};
use std::iter;

const E_BOUND_LEN: usize = (2 * E_BOUND + 1).next_power_of_two().ilog2() as usize;
const K1_BOUND_LEN: usize = (2 * K1_BOUND + 1).next_power_of_two().ilog2() as usize;

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

pub struct BfvEncryptBlock {
    num_reps: usize,
}

impl BfvEncryptBlock {
    pub const LOG2_POLY_SIZE: usize = 10;

    pub const fn log2_size(&self) -> usize {
        // self.perm.log2_size()
        Self::LOG2_POLY_SIZE + 1
    }

    pub const fn r2i_bound_len() -> usize {
        (2 * R2_BOUNDS[0] + 1).next_power_of_two().ilog2() as usize
    }

    pub fn r1i_bound_len(&self) -> Vec<usize> {
        R1_BOUNDS
            .into_iter()
            .take(self.num_reps)
            .map(|b| (2 * b + 1).next_power_of_two().ilog2() as usize)
            // .max()
            .collect()
    }

    // single block
    pub fn configure<F: PrimeField, E: ExtensionField<F>>(
        &self,
        circuit: &mut Circuit<F, E>,
        s: NodeId,
        e: NodeId,
        k1: NodeId,
    ) -> NodeId {
        let poly_log2_size = Self::LOG2_POLY_SIZE;
        let log2_size = self.log2_size();

        let ai = circuit.insert(InputNode::new(log2_size, self.num_reps));

        let r1is = iter::repeat_with(|| circuit.insert(InputNode::new(log2_size, 1)))
            .take(self.num_reps)
            .collect_vec();

        for i in 0..self.num_reps {
            let log2_t_size = self.r1i_bound_len()[i];
            println!("log2_t_size {:?}", log2_t_size);
            let r1i_m = circuit.insert(InputNode::new(log2_t_size, 1));
            let r1i_t = circuit.insert(InputNode::new(log2_t_size, 1));
            let r1i_range = circuit.insert(LogUpNode::new(log2_t_size, log2_size, 1));
            let r1i = r1is[i];

            connect!(circuit {
                r1i_range <- r1i_m, r1i_t, r1i;
            });
        }

        let r1is_par = {
            let r1i_size = 1usize << log2_size;
            println!("r1i_size {:?}", r1i_size);
            // or w/0 cycle since num reps int node??
            let gates = (0..self.num_reps)
                .flat_map(|i| (0..r1i_size).map(move |j| VanillaGate::relay((i, j))))
                .collect_vec();
            println!("r1is_par gates {:?}", gates.len());

            circuit.insert(VanillaNode::new(
                self.num_reps,
                log2_size,
                gates.clone(),
                1,
            ))
        };

        r1is.iter()
            .take(self.num_reps)
            .for_each(|&r1i| circuit.connect(r1i, r1is_par));

        let r2i = circuit.insert(InputNode::new(poly_log2_size, self.num_reps));

        let s_eval = circuit.insert(FftNode::forward(log2_size + self.num_reps - 1));
        let ai_eval = circuit.insert(FftNode::forward(log2_size + self.num_reps - 1));
        let sai_eval = {
            let gates = (0..1usize << log2_size)
                .map(|i| VanillaGate::mul((0, i), (1, i)))
                .collect_vec();
            circuit.insert(VanillaNode::new(2, log2_size, gates, self.num_reps))
        };

        let sai = circuit.insert(FftNode::inverse(log2_size + self.num_reps - 1));

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
                .map(|i| {
                    VanillaGate::new(
                        None,
                        vec![
                            (None, (0, i)),
                            (None, (1, i)),
                            (Some(F::from_str_vartime(K0IS[0]).unwrap()), (2, i)),
                            (Some(F::from_str_vartime(QIS[0]).unwrap()), (3, i)),
                            (None, (4, i)),
                        ],
                        vec![],
                    )
                })
                .collect();
            circuit.insert(VanillaNode::new(5, log2_size, gates, self.num_reps))
        };

        connect!(circuit {
            s_eval <- s;
            ai_eval <- ai;
            sai_eval <- s_eval, ai_eval;
            sai <- sai_eval;
            r2i_cyclo <- r2i;
            sum <- sai, e, k1, r1is_par, r2i_cyclo;
        });

        r1is_par
    }
}

pub struct BfvEncrypt {
    block: BfvEncryptBlock,
}

impl BfvEncrypt {
    pub fn new(num_reps: usize) -> Self {
        Self {
            block: BfvEncryptBlock { num_reps },
        }
    }

    pub const LOG2_POLY_SIZE: usize = 10;

    pub const fn log2_size(&self) -> usize {
        // self.perm.log2_size()
        Self::LOG2_POLY_SIZE + 1
    }

    pub fn configure<F: PrimeField, E: ExtensionField<F>>(
        &self,
        circuit: &mut Circuit<F, E>,
    ) -> NodeId {
        let log2_size = self.log2_size();

        let s = circuit.insert(InputNode::new(log2_size, self.block.num_reps));
        let e = circuit.insert(InputNode::new(log2_size, self.block.num_reps));
        let k1 = circuit.insert(InputNode::new(log2_size, self.block.num_reps));

        let s_m = circuit.insert(InputNode::new(2, 1));
        let s_t = circuit.insert(InputNode::new(2, 1));
        let s_range = circuit.insert(LogUpNode::new(2, log2_size + self.block.num_reps - 1, 1));

        let e_m = circuit.insert(InputNode::new(E_BOUND_LEN, 1));
        let e_t = circuit.insert(InputNode::new(E_BOUND_LEN, 1));
        let e_range = circuit.insert(LogUpNode::new(
            E_BOUND_LEN,
            log2_size + self.block.num_reps - 1,
            1,
        ));

        let k1_m = circuit.insert(InputNode::new(K1_BOUND_LEN, 1));
        let k1_t = circuit.insert(InputNode::new(K1_BOUND_LEN, 1));
        let k1_range = circuit.insert(LogUpNode::new(
            K1_BOUND_LEN,
            log2_size + self.block.num_reps - 1,
            1,
        ));

        connect!(circuit {
            s_range <- s_m, s_t, s;
            e_range <- e_m, e_t, e;
            k1_range <- k1_m, k1_t, k1;
        });

        self.block.configure(circuit, s, e, k1)
    }

    pub fn gen_values<F: PrimeField, E: ExtensionField<F>>(
        &self,
        args: BfvSkEncryptArgs,
    ) -> Vec<BoxMultilinearPoly<'static, F, E>> {
        let log2_size = self.log2_size();
        println!("log2_size {:?}", log2_size);

        let s = Poly::<F>::new_padded(args.s.clone(), log2_size);
        let e = Poly::<F>::new_shifted(args.e.clone(), (1 << log2_size) - 1);
        let k1 = Poly::<F>::new_shifted(args.k1.clone(), (1 << log2_size) - 1);

        let mut r2is = vec![];
        let mut r1is = vec![];
        let mut ais = vec![];
        let mut ct0is = vec![];

        let mut qi_constants = vec![];
        let mut k0i_constants = vec![];

        for z in 0..min(args.ct0is.len(), self.block.num_reps) {
            let r2i = Poly::<F>::new(args.r2is[z].clone());
            r2is.push(r2i.as_ref().to_vec());

            let r1i = Poly::<F>::new_padded(args.r1is[z].clone(), log2_size);
            r1is.push(r1i.as_ref().to_vec());

            let ai = Poly::<F>::new_padded(args.ais[z].clone(), log2_size);
            ais.extend(ai.as_ref().to_vec());

            let ct0i = Poly::<F>::new_shifted(args.ct0is[z].clone(), 1 << log2_size);
            let mut ct0i = ct0i.as_ref()[1..].to_vec();
            ct0i.push(F::ZERO);
            ct0is.extend(ct0i);

            qi_constants.push(F::from_str_vartime(QIS[z]).unwrap());
            k0i_constants.push(F::from_str_vartime(K0IS[z]).unwrap());
        }

        let ss = (0..self.block.num_reps)
            .flat_map(|_| s.as_ref().to_vec())
            .collect_vec();
        let es = (0..self.block.num_reps)
            .flat_map(|_| e.as_ref().to_vec())
            .collect_vec();
        let k1s = (0..self.block.num_reps)
            .flat_map(|_| k1.as_ref().to_vec())
            .collect_vec();

        let omega = root_of_unity(log2_size + self.block.num_reps - 1);

        let s_eval = {
            let mut buf = ss.clone();
            println!("s_eval size {:?}", buf.len());
            radix2_fft(&mut buf, omega);
            buf
        };

        let ai_eval = {
            let mut buf = ais.clone();
            println!("ai_eval size {:?}", buf.len());
            radix2_fft(&mut buf, omega);
            buf
        };

        let sai_eval = izip!(s_eval.iter(), ai_eval.iter())
            .map(|(s, ai)| *s * *ai)
            .collect_vec();

        println!("sai_eval size {:?}", sai_eval.len());

        let sai = {
            let mut buf = sai_eval.clone();
            radix2_ifft(&mut buf, omega);
            buf
        };

        let r2is_cyclo = r2is
            .iter()
            .take(self.block.num_reps)
            .flat_map(|r2i| {
                let mut result = vec![F::ZERO; 2 * N]; // Allocate result vector of size 2N-1

                println!("r2i size {:?}", r2i.len());

                for i in 0..r2i.len() {
                    result[i] += r2i[i]; // Add P(x)
                    result[i + N] += r2i[i]; // Add P(x) * x^N
                }
                result
            })
            .collect_vec();

        let r2is = r2is
            .into_iter()
            .take(self.block.num_reps)
            .flat_map(|mut r2i| {
                r2i.push(F::ZERO);
                r2i
            })
            .collect_vec();

        let ct0i_check = sai
            .iter()
            .zip(es.iter())
            .zip(k1s.iter())
            .zip(r1is.iter().flatten())
            .zip(r2is_cyclo.iter())
            .map(|((((&sai0, &e), &k1), &r1i), &r2i_cyclo)| {
                sai0 + e + k1 * k0i_constants[0] + r1i * qi_constants[0] + r2i_cyclo
            })
            .collect_vec();

        let (s_t, s_m) = {
            let t = [F::ZERO, F::ONE, F::ZERO - F::ONE, F::ZERO];
            let mut m = [F::ZERO; 4];
            ss.iter().for_each(|s| {
                if let Some(i) = t.iter().position(|e| e == s) {
                    m[i] += F::ONE;
                }
            });

            (t.to_vec(), m.to_vec())
        };

        let (e_t, e_m) = {
            let mut t = (0..=E_BOUND)
                .flat_map(|b| [F::ZERO - F::from(b), F::from(b)])
                .collect_vec();
            t.resize(1 << E_BOUND_LEN, F::ZERO);
            let mut m = vec![F::ZERO; 1 << E_BOUND_LEN];
            es.iter().for_each(|s| {
                if let Some(i) = t.iter().position(|e| e == s) {
                    m[i] += F::ONE; // TODO: how do we infore multiplicities in the GKR circuit?
                }
            });

            (t, m)
        };

        let (k1_t, k1_m) = {
            let mut t = (0..=K1_BOUND)
                .flat_map(|b| [F::ZERO - F::from(b), F::from(b)])
                .collect_vec();
            t.resize(1 << K1_BOUND_LEN, F::ZERO);
            let mut m = vec![F::ZERO; 1 << K1_BOUND_LEN];
            k1s.iter().for_each(|s| {
                if let Some(i) = t.iter().position(|e| e == s) {
                    m[i] += F::ONE;
                }
            });

            (t, m)
        };

        // let (r2i_t, r2i_m) = {
        //     let mut t = (0..=max(R2_BOUNDS[0], R2_BOUNDS[1]))
        //         .flat_map(|b| [F::ZERO - F::from(b), F::from(b)])
        //         .collect_vec();
        //     t.resize(1 << R2I_BOUND_LEN, F::ZERO);
        //     let mut m = vec![F::ZERO; 1 << R2I_BOUND_LEN];
        //     r2is.iter().for_each(|s| {
        //         if let Some(i) = t.iter().position(|e| e == s) {
        //             m[i] += F::ONE;
        //         }
        //     });

        //     (t, m)
        // };

        // this is not ideal, need to have a separate range for each parallel sub circuit
        let r1i_range_values = izip!(R1_BOUNDS, self.block.r1i_bound_len(), &r1is)
            .map(|(bound, bound_len, r1i)| {
                let mut t = (0..=bound)
                    .flat_map(|b| [F::ZERO - F::from(b), F::from(b)])
                    .collect_vec();
                t.resize(1 << bound_len, F::ZERO);
                let mut m = vec![F::ZERO; 1 << bound_len];
                r1i.iter().for_each(|s| {
                    if let Some(i) = t.iter().position(|e| e == s) {
                        m[i] += F::ONE;
                    }
                });

                (m, t)
            })
            .take(self.block.num_reps)
            .flat_map(|(m, t)| [m, t, vec![F::ZERO]])
            .collect_vec();

        let r1is_par = r1is
            .iter()
            .take(self.block.num_reps)
            .flatten()
            .cloned()
            .collect_vec();

        println!(
            "r1is sizes {:?}",
            r1is.iter().map(|r1i| r1i.len()).collect_vec()
        );

        println!("r1is_par size {:?}", r1is_par.len());

        chain_par![
            [ss, es, k1s],
            [s_m, s_t, vec![F::ZERO]],   // s_range
            [e_m, e_t, vec![F::ZERO]],   // e_range
            [k1_m, k1_t, vec![F::ZERO]], // k1_range
            [ais],
            r1is,
            r1i_range_values, // r1i_range
            [r1is_par],
            [r2is],
            [s_eval, ai_eval],
            [sai_eval],
            [sai],
            [r2is_cyclo],
            [ct0i_check]
        ]
        .map(box_dense_poly)
        .collect()
    }

    // fn prove<F: PrimeField, E: ExtensionField<F>>(args: BfvSkEncryptArgs, pcs: impl PCS<F, E>,
    //     transcript: &mut impl TranscriptWrite<F, E>,) {
    //     let bfv = BfvEncrypt::new(2);

    //     let circuit = {
    //         let mut circuit = Circuit::default();
    //         bfv.configure(&mut circuit);
    //         circuit
    //     };
    //     let values = bfv.gen_values(args, pcs);

    //     // let values = circuit.evaluate(expected_values);
    //     // assert_polys_eq(&values, &expected_values);
    //     gkr::prove_gkr(&circuit, &values, output_claims, transcript);

    // }
}

pub fn bfv_encrypt_circuit<F: PrimeField + From<u64>, E: ExtensionField<F>>(
    bfv: BfvEncrypt,
    args: BfvSkEncryptArgs,
) -> (Circuit<F, E>, Vec<BoxMultilinearPoly<'static, F, E>>) {
    let circuit = {
        let mut circuit = Circuit::default();
        bfv.configure(&mut circuit);
        circuit
    };
    let values = bfv.gen_values(args);

    (circuit, values)
}

fn root_of_unity<F: PrimeField>(k: usize) -> F {
    assert!(k <= F::S as usize);
    squares(F::ROOT_OF_UNITY).nth(F::S as usize - k).unwrap()
}

pub fn radix2_ifft<F: PrimeField>(buf: &mut [F], omega: F) {
    let n = buf.len();
    let omega_inv = omega.invert().unwrap();

    radix2_fft(buf, omega_inv);

    // Normalize the result by dividing by n
    let n_inv = F::from(n as u64).invert().unwrap(); // Assuming `Field` has an inverse method
    buf.iter_mut().for_each(|x| *x *= &n_inv);
}

#[cfg(test)]
mod test {
    use super::*;
    use gkr::{
        dev::run_gkr_with_values,
        util::dev::seeded_std_rng,
        util::{arithmetic::Field, Itertools},
    };
    use goldilocks::{Goldilocks, GoldilocksExt2};
    use std::{fs::File, io::Read};

    #[test]
    pub fn test_sk_enc_valid() {
        let mut rng = seeded_std_rng();

        // 1. Define the inputs of the circuit
        let file_path = "src/data/sk_enc_1024_2x55_65537.json";
        let mut file = File::open(file_path).unwrap();
        let mut data = String::new();
        file.read_to_string(&mut data).unwrap();
        let bfv = BfvEncrypt::new(2);
        let args = serde_json::from_str::<BfvSkEncryptArgs>(&data).unwrap();

        let (circuit, values) = bfv_encrypt_circuit::<Goldilocks, GoldilocksExt2>(bfv, args);
        // let values = circuit.evaluate(expected_values);
        // assert_polys_eq(&values, &expected_values);
        run_gkr_with_values(&circuit, &values, &mut rng);
    }

    #[test]
    pub fn test_fft_linear() {
        type F = Goldilocks;
        fn pad_to_power_of_two(coeffs: Vec<F>, target_len: usize) -> Vec<F> {
            let mut padded_coeffs = coeffs;
            padded_coeffs.resize(target_len, F::ZERO);
            padded_coeffs
        }

        // Define two pairs of polynomials with 4 coefficients each
        let p1_coeffs: Vec<F> = [1u128, 2, 3, 4]
            .iter()
            .map(|&x| F::from_u128(x))
            .collect_vec();
        let q1_coeffs: Vec<F> = [5u128, 6, 7, 8]
            .iter()
            .map(|&x| F::from_u128(x))
            .collect_vec();
        let p2_coeffs: Vec<F> = [9u128, 10, 11, 12]
            .iter()
            .map(|&x| F::from_u128(x))
            .collect_vec();
        let q2_coeffs: Vec<F> = [13u128, 14, 15, 16]
            .iter()
            .map(|&x| F::from_u128(x))
            .collect_vec();

        // Assume we want to perform FFT over size 16 (next power of two for length 8 + 8)
        let log2_size = 4; // log2(16) = 4
        let omega = root_of_unity(log2_size);

        // Approach 2: concatenate then pad
        // let p_coeffs = [p1_coeffs.clone(), q1_coeffs.clone()].concat();
        // let q_coeffs = [p2_coeffs.clone(), q2_coeffs.clone()].concat();

        let mut p_coeffs = [p1_coeffs.clone(), p2_coeffs.clone()].concat();
        let mut q_coeffs = [q1_coeffs.clone(), q2_coeffs.clone()].concat();

        // Pad polynomials to length 8 (next power of two for 4 + 4 coefficients)
        let padded_p1 = pad_to_power_of_two(p1_coeffs, 8);
        let padded_q1 = pad_to_power_of_two(q1_coeffs, 8);
        let padded_p2 = pad_to_power_of_two(p2_coeffs, 8);
        let padded_q2 = pad_to_power_of_two(q2_coeffs, 8);

        // bit_reverse_permute(&mut p_coeffs);
        // bit_reverse_permute(&mut q_coeffs);

        // Approach 1: Concatenate padded coefficients for batched FFT
        // let p_coeffs = [padded_p1.clone(), padded_p2.clone()].concat();
        // let q_coeffs = [padded_q1.clone(), padded_q2.clone()].concat();

        // let p_coeffs = [padded_p1.clone(), padded_q1.clone()].concat();
        // let q_coeffs = [padded_p2.clone(), padded_q2.clone()].concat();

        // Approach 2: Concatenate padded coefficients for batched FFT
        let p_coeffs = pad_to_power_of_two(p_coeffs, 16);
        let q_coeffs = pad_to_power_of_two(q_coeffs, 16);

        // Perform FFT on the concatenated coefficients
        let mut fft_p = p_coeffs.clone();
        let mut fft_q = q_coeffs.clone();
        radix2_fft(&mut fft_p, omega);
        radix2_fft(&mut fft_q, omega);

        println!("fft_p {:?}", fft_p);
        println!("fft_q {:?}", fft_q);

        // Perform element-wise multiplication in the frequency domain
        let fft_pq: Vec<F> = fft_p
            .iter()
            .zip(fft_q.iter())
            .map(|(p, q)| *p * *q)
            .collect();

        // Perform IFFT to get the resulting coefficients in the time domain
        let mut pq_coeffs = fft_pq.clone();
        radix2_ifft(&mut pq_coeffs, omega);

        // Split the result back into two parts
        let pq1_coeffs = &pq_coeffs[..8]; // Result for P1 * Q1
        let pq2_coeffs = &pq_coeffs[8..]; // Result for P2 * Q2

        // To verify, perform individual FFT-based multiplications for P1 * Q1 and P2 * Q2
        let omega = root_of_unity(log2_size - 1);

        let pq1_check = {
            let mut fft_p1 = padded_p1.clone();
            let mut fft_q1 = padded_q1.clone();
            radix2_fft(&mut fft_p1, omega);
            radix2_fft(&mut fft_q1, omega);

            let fft_pq1: Vec<F> = fft_p1
                .iter()
                .zip(fft_q1.iter())
                .map(|(p, q)| *p * *q)
                .collect();
            let mut pq1_check = fft_pq1.clone();
            radix2_ifft(&mut pq1_check, omega);
            pq1_check
        };

        let pq2_check = {
            let mut fft_p2 = padded_p2.clone();
            let mut fft_q2 = padded_q2.clone();
            radix2_fft(&mut fft_p2, omega);
            radix2_fft(&mut fft_q2, omega);

            // println!("fft_p2 {:?}", fft_p2);
            // println!("fft_q2 {:?}", fft_q2);

            let fft_pq2: Vec<F> = fft_p2
                .iter()
                .zip(fft_p2.iter())
                .map(|(p, q)| *p * *q)
                .collect();
            let mut pq2_check = fft_pq2.clone();
            radix2_ifft(&mut pq2_check, omega);
            pq2_check
        };

        let pq_check = [pq1_check, pq2_check].concat();

        // Check that the batched result matches the individual results
        // assert_eq!(pq1_coeffs, pq1_check, "PQ1 does not match");
        // assert_eq!(pq2_coeffs, pq2_check, "PQ2 does not match");

        assert_eq!(pq_coeffs, pq_check, "PQ does not match");
    }

    fn bit_reverse_permute<F: Clone>(coeffs: &mut [F]) {
        let n = coeffs.len();
        let mut i = 0;
        for j in 1..n {
            let mut bit = n >> 1;
            while i >= bit {
                i -= bit;
                bit >>= 1;
            }
            i += bit;
            if j < i {
                coeffs.swap(j, i);
            }
        }
    }
}
