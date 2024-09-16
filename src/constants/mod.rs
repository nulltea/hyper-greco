mod sk_enc_constants_1024_1x27_65537;
mod sk_enc_constants_16384_8x54_65537;
mod sk_enc_constants_2048_1x52_65537;
mod sk_enc_constants_32768_15x59_65537;
mod sk_enc_constants_4096_2x55_65537;
mod sk_enc_constants_8192_4x55_65537;

pub use sk_enc_constants_1024_1x27_65537::SkEnc1024_1x27_65537;
pub use sk_enc_constants_16384_8x54_65537::SkEnc16384_8x54_65537;
pub use sk_enc_constants_2048_1x52_65537::SkEnc2048_1x52_65537;
pub use sk_enc_constants_32768_15x59_65537::SkEnc32768_15x59_65537;
pub use sk_enc_constants_4096_2x55_65537::SkEnc4096_2x55_65537;
pub use sk_enc_constants_8192_4x55_65537::SkEnc8192_4x55_65537;

pub trait BfvSkEncryptConstans<const K: usize> {
    /// `N` is the degree of the cyclotomic polynomial defining the ring `Rq = Zq[X]/(X^N + 1)`.
    const N: usize;
    /// `N_LOG2` is the logarithm of `N`.
    const N_LOG2: usize = Self::N.ilog2() as usize;
    /// The coefficients of the polynomial `e` should exist in the interval `[-E_BOUND, E_BOUND]` where `E_BOUND` is the upper bound of the gaussian distribution with 𝜎 = 3.2
    const E_BOUND: u64;
    /// The coefficients of `k1` should exist in the interval `[-K1_BOUND, K1_BOUND]` where `K1_BOUND` is equal to `(t-1)/2`
    const K1_BOUND: u64;
    /// The coefficients of the polynomial `s` should exist in the interval `[-S_BOUND, S_BOUND]`.
    const S_BOUND: u64;
    /// The coefficients of the polynomials `r1is` should exist in the interval `[-R1_BOUND[i], R1_BOUND[i]]` where `R1_BOUND[i]` is equal to `(qi-1)/2`
    const R1_BOUNDS: [u64; K];
    /// The coefficients of the polynomials `r2is` should exist in the interval `[-R2_BOUND[i], R2_BOUND[i]]` where `R2_BOUND[i]` is equal to $\frac{(N+2) \cdot \frac{q_i - 1}{2} + B + \frac{t - 1}{2} \cdot |K_{0,i}|}{q_i}$
    const R2_BOUNDS: [u64; K];
    /// List of scalars `qis` such that `qis[i]` is the modulus of the i-th CRT basis of `q` (ciphertext space modulus)
    const QIS: [&str; K];
    /// List of scalars `k0is` such that `k0i[i]` is equal to the negative of the multiplicative inverses of t mod qi.
    const K0IS: [&str; K];
}
