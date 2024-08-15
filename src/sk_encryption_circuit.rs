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
use rayon::iter::ParallelIterator;
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

pub struct BfvBlockInput {
    s: NodeId,
    ai: NodeId,
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

    // pub fn alloc_input<F: PrimeField, E: ExtensionField<F>>(
    //     &self,
    //     circuit: &mut Circuit<F, E>,
    // ) -> BfvBlockInput {
    //     let s = circuit.insert(InputNode::new(Self::LOG2_POLY_SIZE, 1));
    //     let ai = circuit.insert(InputNode::new(Self::LOG2_POLY_SIZE, 1));
    //     BfvBlockInput { s, ai }
    // }

    pub fn configure<F: PrimeField, E: ExtensionField<F>>(&self, circuit: &mut Circuit<F, E>) {
        // single block

        let log2_size = self.log2_size();

        let s = circuit.insert(InputNode::new(log2_size, 1));
        let ai = circuit.insert(InputNode::new(log2_size, 1));
        let e = circuit.insert(InputNode::new(log2_size, 1));
        let k1 = circuit.insert(InputNode::new(log2_size, 1));
        let r1i = circuit.insert(InputNode::new(log2_size, 1));
        let r2i_cyclo = circuit.insert(InputNode::new(log2_size, 1));
        // let cyclo = circuit.insert(InputNode::new(log2_size, 1));

        let s_eval = circuit.insert(FftNode::forward(log2_size));
        let ai_eval = circuit.insert(FftNode::forward(log2_size));
        let mul = {
            let gates = vec![VanillaGate::mul((0, 0), (1, 0))];
            circuit.insert(VanillaNode::new(2, 0, gates, 1 << log2_size))
        }; // why this works with log2_sub_input_size = 0 and num reps = 2^log2_input_size

        let sai0 = circuit.insert(FftNode::inverse(log2_size));

        // let mul_cyclo = {
        //     let gates = vec![VanillaGate::mul((0, 0), (1, 0))];
        //     circuit.insert(VanillaNode::new(2, 0, gates, 1 << log2_size))
        // };

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

        // let sum = {
        //     let gates = vec![
        //         VanillaGate::new(None, vec![(None, (0, 0)), (None, (1, 0)), (Some(F::from_str_vartime(K0IS[0]).unwrap()), (2, 0)),  /*(Some(F::from_str_vartime(QIS[0]).unwrap()), (3, 0))*/], vec![]),
        //     ];
        //     circuit.insert(VanillaNode::new(3, 0, gates, 1 << log2_size))
        // };

        // let sum = {
        //     let gates = vec![VanillaGate::sum(vec![(0, 0), (1, 0)])];
        //     circuit.insert(VanillaNode::new(2, 1, gates, 1 << log2_size))
        // };

        connect!(circuit {
            s_eval <- s;
            ai_eval <- ai;
            mul <- s_eval, ai_eval;
            sai0 <- mul;
            sum <- sai0, e, k1, r1i, r2i_cyclo;
        });
    }

    pub fn gen_values<F: PrimeField, E: ExtensionField<F>>(
        &self,
        args: BfvSkEncryptArgs,
    ) -> Vec<BoxMultilinearPoly<'static, F, E>> {
        let log2_size = self.log2_size();
        println!("log2_size {:?}", log2_size);
        // let rate = bfv.rate();

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
            ct0is.push(ct0i);
        }

        let omega = root_of_unity(log2_size);

        // println!("s {:?}", &s.coefficients[..10]);
        // println!("e {:?}", &e.coefficients[..10]);

        let s_eval = {
            let mut buf = s.as_ref().to_vec();
            radix2_fft(&mut buf, omega);
            buf
        };

        let ai0_eval = {
            let mut buf = ais[0].as_ref().to_vec();
            // buf.resize(1 << (log2_size + 1), F::ZERO);
            radix2_fft(&mut buf, omega);
            buf
        };

        // println!("-------------------");
        // // println!("s_eval[ :10] {:?}", &s_eval[..10]);
        // println!("s_eval[-10:] {:?}", &s_eval[s_eval.len() - 11..]);
        // // println!("ai0_eval[ :10] {:?}", &ai0_eval[..10]);
        // println!("ai0_eval[-10:] {:?}", &ai0_eval[ai0_eval.len() - 11..]);

        // println!("-------------------");

        let sai0 = s_eval
            .iter()
            .zip(ai0_eval.iter())
            .map(|(s, ai)| *s * *ai)
            .collect_vec();

        let mut sai0_poly = {
            let mut buf = sai0.clone();
            radix2_ifft(&mut buf, omega);
            buf
        };

        // sai0_poly.pop();
        // sai0_poly.insert(0, F::ZERO);

        println!("sai0[ :10] {:?}", &sai0_poly[..10]);
        println!("sai0[-10:] {:?}", &sai0_poly[sai0_poly.len()-11..]);
        println!("-------------------");

        let cyclo = Poly::<F>::cyclo_padded(log2_size);

        // let r2i_cyclo = {
        //     let mut coeffs = r2is[0].as_ref().to_vec();
        //     println!("r2i len {:?}", args.r2is[0].len());
        //     let lbs = coeffs[N-1..].iter().cloned().collect_vec();
        //     coeffs.splice(N-1.., lbs);
        //     coeffs
        // };

        let r2i_cyclo = {
            // println!("r2i len {:?}", args.r2is[0].len());
            // println!("r2i*c len {:?}", N + N - 1);

            let r2i0 = Poly::<F>::new(args.r2is[0].clone());
            let mut result = vec![F::ZERO; N + N - 1]; // Allocate result vector of size 2N-1

            // Step 1: Add P(x) to result
            for i in 0..args.r2is[0].len() {
                result[i] += r2i0.coefficients[i]; // Add P(x)
                result[i + N] += r2i0.coefficients[i]; // Add P(x) * x^N
            }
            result.push(F::ZERO);
            result
        };

        // println!("r2i_cyclo[ :10] {:?}", &r2i_cyclo[..10]);
        // println!("r2i_cyclo[-10:] {:?}", &r2i_cyclo[N-2..N+2]);

        // {
        //     let mut coeffs = vec![0, 0, 0, 0, 0, 1, 2, 3];
        //     let lbs = coeffs[..4].iter().cloned().collect_vec();
        //     coeffs.splice(4-1.., lbs);
        //     println!("coeffs test {:?}", coeffs);
        // };

        // println!("r2i_cyclo coeefffs {:?}", &r2i_cyclo[r2i_cyclo.len() - 50..]);

      


        let r1iqis = r1is[0]
            .as_ref()
            .iter()
            .map(|r1i| *r1i * qi_constants[0])
            .collect_vec();


        // *sai0 + *e + *k1 * k0i_constants[0] + *r1i * qi_constants[0] + *r2i_cyclo
        // let k1k0i = k1
        //     .as_ref()
        //     .iter()
        //     .map(|k1| *k1 * k0i_constants[0])
        //     .collect_vec();

        // println!("k1k0i[ :10] {:?}", &k1k0i[..10]);
        // println!("k1k0i[-10:] {:?}", &k1k0i[k1k0i.len() - 10..]);

        let ct_hat = sai0_poly
        .iter()
        .zip(e.as_ref().iter())
        .zip(k1.as_ref().iter())
        .zip(r1is[0].as_ref().iter())
        .zip(r2i_cyclo.iter())
        .map(|((((sai0, e), k1), r1i), r2i_cyclo)| {
            *sai0 + *e
        })
        .collect_vec();

        // println!("ct_hat[ :10] {:?}", &ct_hat[..10]);
        println!("e len {:?}", e.coefficients.len());
        println!("e[-10:] {:?}", &e.coefficients[e.coefficients.len()-11..]);
        println!("ct_hat[ :10] {:?}", &ct_hat[..10]);
        println!("ct_hat[-10:] {:?}", &ct_hat[ct_hat.len()-11..]);
        println!("-------------------");


        println!("r1i[ :10] {:?}", &r1is[0].coefficients[..10]);
        println!("r1iqis[ :10] {:?}", &r1iqis[..10]);
        println!("-------------------");
        // println!("r1iqis[-10:] {:?}", &r1iqis[r1iqis.len() - 10..]);

        let sai0_e = sai0_poly
            .iter()
            .zip(e.as_ref().iter())
            .zip(k1.as_ref().iter())
            .zip(r1is[0].as_ref().iter())
            .zip(r2i_cyclo.iter())
            .map(|((((sai0, e), k1), r1i), r2i_cyclo)| {
                *sai0 + *e + *k1 * k0i_constants[0] + *r1i * qi_constants[0] + *r2i_cyclo
            })
            .collect_vec();
        //*sai0 + *e + *k1 * k0i_constants[0] + *r1i * qi_constants[0] + *r2i_cyclo

        println!("rhs[ :10] {:?}", &sai0_e[..10]);
        println!("rhs[-10:] {:?}", &sai0_e[sai0_e.len() - 11..]);
        println!("-------------------");

        // let k1_k01 = k1.as_ref().iter().map(|k1| *k1 * k01).collect_vec();

        // println!("r1is 0 {:?}", r1is[0].coefficients.len());
        // println!("r2i_cyclo {:?}", cyclo.coefficients.len());

        
        let mut ct0i = ct0is[0].as_ref()[1..].to_vec();
        println!("ct0i len {:?}", ct0i.len());
        ct0i.push(F::ZERO);

        println!("cti[ :10] {:?}", &ct0i[..10]);
        println!("cti[-10:] {:?}", &ct0i[ct0i.len() - 11..]);

        // println!("sai0_e coeefffs {:?}", &sai0_e[sai0_e.len() - 100..]);

        chain_par![
            [
                s,
                ais[0].clone(),
                e,
                k1,
                r1is[0].clone(),
            ]
            .map(|p| p.as_ref().to_vec()),
            [r2i_cyclo],
            [s_eval, ai0_eval],
            [sai0],
            [sai0_poly],
            // [k1_k01],
            [ct0i],
        ]
        .map(box_dense_poly)
        .collect()
    }
}

pub fn bfv_encrypt_circuit<F: PrimeField + From<u64>, E: ExtensionField<F>>(
    bfv: BfvEncryptBlock,
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
        util::dev::{assert_polys_eq, seeded_std_rng},
    };
    use goldilocks::{Goldilocks, GoldilocksExt2};
    use halo2_curves::bn256;
    use std::{fs::File, io::Read};

    #[test]
    pub fn test_sk_enc_valid() {
        let mut rng = seeded_std_rng();

        // 1. Define the inputs of the circuit
        let file_path = "src/data/sk_enc_1024_2x55_65537.json";
        let mut file = File::open(file_path).unwrap();
        let mut data = String::new();
        file.read_to_string(&mut data).unwrap();
        let bfv = BfvEncryptBlock { num_reps: 1 };
        let args = serde_json::from_str::<BfvSkEncryptArgs>(&data).unwrap();

        let (circuit, values) = bfv_encrypt_circuit::<bn256::Fr, bn256::Fr>(bfv, args);
        // let values = circuit.evaluate(expected_values);
        // assert_polys_eq(&values, &expected_values);
        run_gkr_with_values(&circuit, &values, &mut rng);
    }
}
