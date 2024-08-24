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

// pub struct PolyAssigned<F: ScalarField> {
//     pub assigned_coefficients: Vec<AssignedValue<F>>,
// }

// impl<F: ScalarField> PolyAssigned<F> {
//     pub fn new(ctx: &mut Context<F>, poly: Poly<F>) -> Self {
//         let assigned_coefficients = ctx.assign_witnesses(poly.coefficients);
//         PolyAssigned {
//             assigned_coefficients,
//         }
//     }

//     /// Adds `upper_bound` to the coefficients of the polynomial and constrains them to be in the range `[0, 2*upper_bound]`.
//     pub fn range_check(&self, ctx_gate: &mut Context<F>, range: &RangeChip<F>, upper_bound: u64) {
//         let bound_constant = Constant(F::from(upper_bound));

//         for coeff in &self.assigned_coefficients {
//             let shifted_coeff = range.gate().add(ctx_gate, *coeff, bound_constant);
//             range.check_less_than_safe(ctx_gate, shifted_coeff, (2 * upper_bound) + 1);
//         }
//     }

//     pub fn enforce_eval_at_gamma(
//         &self,
//         ctx_rlc: &mut Context<F>,
//         rlc: &RlcChip<F>,
//     ) -> AssignedValue<F> {
//         let rlc_trace = rlc.compute_rlc_fixed_len(ctx_rlc, self.assigned_coefficients.clone());
//         rlc_trace.rlc_val
//     }
// }
