use gkr::ff_ext::ff::PrimeField;
use itertools::Itertools;

/// Struct to store the coefficients of a polynomial as Field Elements
/// The coefficients are stored starting from the highest degree term
#[derive(Clone, Debug)]
pub struct Poly<F: PrimeField> {
    pub coefficients: Vec<F>,
}

impl<F: PrimeField> Poly<F> {
    pub fn new(coefficients: Vec<String>) -> Self {
        let coefficients = coefficients
            .iter()
            .map(|coeff| F::from_str_vartime(coeff).unwrap())
            .collect();
        Poly { coefficients }
    }

    pub fn new_padded(coefficients: Vec<String>, log2_size: usize) -> Self {
        let mut coefficients = coefficients
            .iter()
            .map(|coeff| F::from_str_vartime(coeff).unwrap())
            .collect_vec();
        coefficients.resize(1 << log2_size, F::ZERO);

        Poly { coefficients }
    }

    pub fn new_shifted(coefficients: Vec<String>, size: usize) -> Self {
        let coefficients = coefficients
            .iter()
            .map(|coeff| F::from_str_vartime(coeff).unwrap())
            .collect_vec();
        let padding_size = size.saturating_sub(coefficients.len());
        
        let mut shifted = vec![F::ZERO; padding_size];
        shifted.extend(coefficients);
        shifted.resize(size.next_power_of_two(), F::ZERO);

        Poly { coefficients: shifted }
    }

    pub fn to_vec(&self) -> Vec<F> {
        self.coefficients.clone()
    }
}

impl<F: PrimeField> AsRef<[F]> for Poly<F> {
    fn as_ref(&self) -> &[F] {
       &self.coefficients
    }
}
