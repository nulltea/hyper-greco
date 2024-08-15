use crate::constants::sk_enc_constants_1024_2x55_65537::{
    E_BOUND, K0IS, K1_BOUND, N, QIS, R1_BOUNDS, R2_BOUNDS, S_BOUND,
};
use crate::poly::Poly;
use gkr::util::arithmetic::radix2_fft;
use gkr::{
    chain_par,
    circuit::{
        connect,
        node::{FftNode, InputNode, LogUpNode, Node, VanillaGate, VanillaNode},
        Circuit, NodeId,
    },
    ff_ext::ff::PrimeField,
    poly::{box_dense_poly, BoxMultilinearPoly},
    util::{
        arithmetic::{powers, squares, ExtensionField, Field},
        chain, izip, Itertools,
    },
};
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};
use serde::Deserialize;

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
    // pub const fn num_bits(&self) -> usize {
    //     self.num_bits
    // }

    // pub const fn rate(&self) -> usize {
    //     self.rate
    // }

    // pub const fn num_reps(&self) -> usize {
    //     self.perm.num_reps()
    // }

    // pub const fn log2_reps(&self) -> usize {
    //     self.perm.log2_reps()
    // }

    pub const fn log2_size(&self) -> usize {
        // self.perm.log2_size()
        Self::LOG2_POLY_SIZE + 1
    }

    // pub fn alloc_input<F: PrimeField, E: ExtensionField<F>>(
    //     &self,
    //     circuit: &mut Circuit<F, E>,
    // ) -> BfvBlockInput {
    //     let s = circuit.insert(InputNode::new(Self::LOG2_POLY_SIZE, 1));
    //     let ai = circuit.insert(InputNode::new(Self::LOG2_POLY_SIZE, 1));
    //     BfvBlockInput { s, ai }
    // }

    pub fn configure<F: PrimeField, E: ExtensionField<F>>(
        &self,
        circuit: &mut Circuit<F, E>,
    ) -> NodeId {
        let log2_size = self.log2_size();

        let s = circuit.insert(InputNode::new(log2_size, 1));
        let e = circuit.insert(InputNode::new(log2_size, 1));
        let k1 = circuit.insert(InputNode::new(log2_size, 1));

        // connect!(circuit {
        //     s_eval <- s;
        //     ai_eval <- ai;
        //     mul <- s_eval, ai_eval;
        //     sai0 <- mul;
        //     sum <- sai0, e, k1, r1i, r2i_cyclo;
        // });

        self.block.configure(circuit, s, e, k1)
    }

    pub fn gen_values<F: PrimeField, E: ExtensionField<F>>(
        &self,
        args: BfvSkEncryptArgs,
    ) -> Vec<BoxMultilinearPoly<'static, F, E>> {
        let log2_size = self.log2_size();
        println!("log2_size {:?}", log2_size);

        let s = Poly::<F>::new_shifted(args.s.clone(), 1 << log2_size);
        let e = Poly::<F>::new_shifted(args.e.clone(), (1 << log2_size) - 1);
        let k1 = Poly::<F>::new_shifted(args.k1.clone(), (1 << log2_size) - 1);

        let mut r2is = vec![];
        let mut r1is = vec![];
        let mut ais = vec![];
        let mut ct0is = vec![];

        let mut qi_constants = vec![];
        let mut k0i_constants = vec![];

        for z in 0..args.ct0is.len() {
            qi_constants.push(F::from_str_vartime(QIS[z]).unwrap());
            k0i_constants.push(F::from_str_vartime(K0IS[z]).unwrap());
        }

        for z in 0..args.ct0is.len() {
            let r2i = Poly::<F>::new_padded(args.r2is[z].clone(), log2_size);
            r2is.push(r2i);

            let r1i = Poly::<F>::new_padded(args.r1is[z].clone(), log2_size);
            r1is.push(r1i);

            let ai = Poly::<F>::new_shifted(args.ais[z].clone(), 1 << log2_size);
            ais.push(ai);

            let ct0i = Poly::<F>::new_shifted(args.ct0is[z].clone(), 1 << log2_size);
            let mut ct0i = ct0i.as_ref()[1..].to_vec();
            ct0i.push(F::ZERO);
            ct0is.push(ct0i);
        }

        let omega = root_of_unity(log2_size);

        let s_eval = {
            let mut buf = s.as_ref().to_vec();
            radix2_fft(&mut buf, omega);
            buf
        };

        let ai0_eval = {
            let mut buf = ais[0].as_ref().to_vec();
            radix2_fft(&mut buf, omega);
            buf
        };

        let sai0 = s_eval
            .iter()
            .zip(ai0_eval.iter())
            .map(|(s, ai)| *s * *ai)
            .collect_vec();

        let sai0_poly = {
            let mut buf = sai0.clone();
            radix2_ifft(&mut buf, omega);
            buf
        };

        let r2i_cyclo = {
            let r2i0 = Poly::<F>::new(args.r2is[0].clone());
            let mut result = vec![F::ZERO; N + N - 1]; // Allocate result vector of size 2N-1

            for i in 0..args.r2is[0].len() {
                result[i] += r2i0.coefficients[i]; // Add P(x)
                result[i + N] += r2i0.coefficients[i]; // Add P(x) * x^N
            }
            result.push(F::ZERO);
            result
        };

        // let ct0i_check = sai0_poly
        //     .iter()
        //     .zip(e.as_ref().iter())
        //     .zip(k1.as_ref().iter())
        //     .zip(r1is[0].as_ref().iter())
        //     .zip(r2i_cyclo.iter())
        //     .map(|((((sai0, e), k1), r1i), r2i_cyclo)| {
        //         *sai0 + *e + *k1 * k0i_constants[0] + *r1i * qi_constants[0] + *r2i_cyclo
        //     })
        //     .collect_vec();

        chain_par![
            [s, e, k1].map(|p| p.as_ref().to_vec()),
            ais.par_iter().take(self.block.num_reps).map(|ai| ai.as_ref().to_vec()),
            r1is.par_iter().take(self.block.num_reps).map(|ai| ai.as_ref().to_vec()),
            [r2i_cyclo],
            [s_eval, ai0_eval],
            [sai0],
            [sai0_poly],
            ct0is.into_par_iter().take(self.block.num_reps),
        ]
        .map(box_dense_poly)
        .collect()
    }
}
pub struct BfvEncryptBlock {
    num_reps: usize,
}

impl BfvEncryptBlock {
    pub const LOG2_POLY_SIZE: usize = 10;
    // pub const fn num_bits(&self) -> usize {
    //     self.num_bits
    // }

    // pub const fn rate(&self) -> usize {
    //     self.rate
    // }

    // pub const fn num_reps(&self) -> usize {
    //     self.perm.num_reps()
    // }

    // pub const fn log2_reps(&self) -> usize {
    //     self.perm.log2_reps()
    // }

    pub const fn log2_size(&self) -> usize {
        // self.perm.log2_size()
        Self::LOG2_POLY_SIZE + 1
    }

    // single block
    pub fn configure<F: PrimeField, E: ExtensionField<F>>(
        &self,
        circuit: &mut Circuit<F, E>,
        s: NodeId,
        e: NodeId,
        k1: NodeId,
    ) -> NodeId {
        let log2_size = self.log2_size();

        let ai = circuit.insert(InputNode::new(log2_size, self.num_reps));
        let r1i = circuit.insert(InputNode::new(log2_size, self.num_reps));
        let r2i_cyclo = circuit.insert(InputNode::new(log2_size, self.num_reps));

        let s_eval = circuit.insert(FftNode::forward(log2_size));
        let ai_eval = circuit.insert(FftNode::forward(log2_size));
        let mul = {
            let gates = vec![VanillaGate::mul((0, 0), (1, 0))];
            circuit.insert(VanillaNode::new(2, 0, gates, 1 << log2_size))
        }; // why this works with log2_sub_input_size = 0 and num reps = 2^log2_input_size

        let sai0 = circuit.insert(FftNode::inverse(log2_size));

        let sum = {
            let gates = vec![VanillaGate::new(
                None,
                vec![
                    (None, (0, 0)),
                    (None, (1, 0)),
                    (Some(F::from_str_vartime(K0IS[0]).unwrap()), (2, 0)),
                    (Some(F::from_str_vartime(QIS[0]).unwrap()), (3, 0)),
                    (None, (4, 0)),
                ],
                vec![],
            )];
            circuit.insert(VanillaNode::new(5, 0, gates, 1 << log2_size))
        };

        connect!(circuit {
            s_eval <- s;
            ai_eval <- ai;
            mul <- s_eval, ai_eval;
            sai0 <- mul;
            sum <- sai0, e, k1, r1i, r2i_cyclo;
        });

        sum
    }
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
    use gkr::{dev::run_gkr_with_values, util::dev::seeded_std_rng};
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
}
