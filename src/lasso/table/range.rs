use gkr::{
    ff_ext::{ff::PrimeField, ExtensionField},
    poly::{MultilinearPolyTerms, PolyExpr},
    util::{arithmetic::inner_product, expression::Expression},
};
use itertools::Itertools;
use std::{iter, marker::PhantomData};

use super::{LassoSubtable, LookupType, SubtableIndices};

#[derive(Clone, Debug, Default)]
pub struct FullLimbSubtable<F, E>(PhantomData<(F, E)>);

impl<F: PrimeField, E: ExtensionField<F>> LassoSubtable<F, E> for FullLimbSubtable<F, E> {
    fn materialize(&self, M: usize) -> Vec<F> {
        (0..M).map(|x| F::from(x as u64)).collect_vec()
    }

    fn evaluate_mle(&self, point: &[E], _: usize) -> E {
        let b = point.len();
        let mut result = E::ZERO;
        for i in 0..b {
            result += point[i] * F::from(1u64 << (i));
        }
        result
    }

    fn evaluate_mle_expr(&self, log2_M: usize) -> MultilinearPolyTerms<F> {
        let limb_init = PolyExpr::Var(0);
        let mut limb_terms = vec![limb_init];
        (1..log2_M).for_each(|i| {
            let coeff = PolyExpr::Pow(Box::new(PolyExpr::Const(F::from(2))), i as u32);
            let x = PolyExpr::Var(i);
            let term = PolyExpr::Prod(vec![coeff, x]);
            limb_terms.push(term);
        });
        MultilinearPolyTerms::new(log2_M, PolyExpr::Sum(limb_terms))
    }
}

impl<F, E> FullLimbSubtable<F, E> {
    pub fn new() -> Self {
        Self(PhantomData)
    }
}

#[derive(Clone, Debug, Default)]
pub struct BoundSubtable<F, E, const BOUND: u64>(PhantomData<(F, E)>);

impl<F: PrimeField, E: ExtensionField<F>, const BOUND: u64> LassoSubtable<F, E>
    for BoundSubtable<F, E, BOUND>
{
    fn materialize(&self, M: usize) -> Vec<F> {
        let bound_bits = BOUND.ilog2() as usize;
        let reminder = 1 << (bound_bits % M.ilog2() as usize);
        let cutoff = (reminder + BOUND % M as u64) as usize;

        (0..M)
            .map(|i| {
                if i < cutoff {
                    F::from(i as u64)
                } else {
                    F::ZERO
                }
            })
            .collect()
    }

    fn evaluate_mle(&self, point: &[E], M: usize) -> E {
        let log2_M = M.ilog2() as usize;
        let b = point.len();

        let bound_bits = BOUND.ilog2() as usize;
        let reminder = 1 << (bound_bits % log2_M);
        let cutoff = reminder + BOUND % (1 << log2_M);
        let cutoff_log2 = cutoff.ilog2() as usize;

        let g_base = 1 << cutoff_log2;
        let num_extra = cutoff - g_base;

        let mut result = E::ZERO;
        for i in 0..b {
            if i < cutoff_log2 {
                result += point[i] * F::from(1u64 << (i));
            } else {
                let mut g_value = E::ZERO;

                if i == cutoff_log2 {
                    for k in 0..num_extra {
                        let mut term = E::from_bases(&[F::from((g_base + k) as u64)]);
                        for j in 0..cutoff_log2 {
                            if (k & (1 << j)) != 0 {
                                term *= point[j];
                            } else {
                                term *= E::ONE - point[j];
                            }
                        }
                        g_value += term;
                    }
                }

                result = (E::ONE - point[i]) * result + point[i] * g_value
            }
        }

        result
    }

    fn evaluate_mle_expr(&self, log2_M: usize) -> MultilinearPolyTerms<F> {
        let M = 1 << log2_M;
        let bound_bits = BOUND.ilog2() as usize;
        let reminder = 1 << (bound_bits % log2_M);
        let cutoff = reminder + BOUND % M;
        let cutoff_log2 = cutoff.ilog2() as usize;

        let rem_init = PolyExpr::Var(0);
        let mut terms = vec![rem_init];
        (1..cutoff_log2).for_each(|i| {
            let coeff = PolyExpr::Pow(Box::new(PolyExpr::Const(F::from(2))), i as u32);
            let x = PolyExpr::Var(i);
            let term = PolyExpr::Prod(vec![coeff, x]);
            terms.push(term);
        });

        let mut result = PolyExpr::Sum(terms);

        let g_base = 1 << cutoff_log2;
        let num_extra = cutoff - g_base;

        (cutoff_log2..log2_M).for_each(|i| {
            if num_extra > 0 && i == cutoff_log2 {
                let mut g_value = PolyExpr::ZERO;
                for k in 0..num_extra {
                    let mut term = PolyExpr::u64((g_base + k) as u64);
                    for j in 0..cutoff_log2 {
                        if (k & (1 << j)) != 0 {
                            term = PolyExpr::mul(term, PolyExpr::Var(j));
                        } else {
                            let t = PolyExpr::sub(PolyExpr::Const(F::ONE), PolyExpr::Var(j));
                            term = PolyExpr::mul(term, t);
                        }
                    }
                    g_value = PolyExpr::add(term, g_value);
                }
                let x = PolyExpr::Var(i);
                let left = PolyExpr::mul(PolyExpr::sub(PolyExpr::ONE, x.clone()), result.clone());
                let right = PolyExpr::mul(x, g_value);
                result = PolyExpr::add(left, right);
            } else {
                let term = PolyExpr::sub(PolyExpr::ONE, PolyExpr::Var(i));
                result = PolyExpr::mul(result.clone(), term);
            }
        });

        MultilinearPolyTerms::new(log2_M, result)
    }
}

impl<F, E, const BOUND: u64> BoundSubtable<F, E, BOUND> {
    pub fn new() -> Self {
        Self(PhantomData)
    }
}

#[derive(Clone, Debug, Default, Copy)]
pub struct RangeLookup<const BOUND: u64>;

impl<const BOUND: u64> LookupType for RangeLookup<BOUND> {
    fn combine_lookups<F: PrimeField>(&self, operands: &[F], _: usize, M: usize) -> F {
        let weight = F::from(M as u64);
        inner_product(
            operands,
            iter::successors(Some(F::ONE), |power_of_weight| {
                Some(*power_of_weight * weight)
            })
            .take(operands.len())
            .collect_vec()
            .iter(),
        )
    }

    fn combine_lookup_expressions<F: PrimeField, E: ExtensionField<F>>(
        &self,
        expressions: Vec<Expression<E, usize>>,
        C: usize,
        M: usize,
    ) -> Expression<E, usize> {
        Expression::distribute_powers(expressions, E::from_bases(&[F::from(M as u64)]))
    }

    // SubtableIndices map subtable to memories
    fn subtables<F: PrimeField, E: ExtensionField<F>>(
        &self,
        C: usize,
        M: usize,
    ) -> Vec<(Box<dyn LassoSubtable<F, E>>, SubtableIndices)> {
        let full = Box::new(FullLimbSubtable::<F, E>::new());
        let log_M = M.ilog2() as usize;
        let bound_bits = BOUND.ilog2() as usize;
        let num_chunks = bound_bits / log_M;
        if bound_bits % log_M == 0 {
            vec![(full, SubtableIndices::from(0..num_chunks))]
        } else {
            let rem = Box::new(BoundSubtable::<F, E, BOUND>::new());
            vec![
                (full, SubtableIndices::from(0..num_chunks)),
                (rem, SubtableIndices::from(num_chunks)),
            ]
        }
    }

    fn output<F: PrimeField>(&self, index: &F) -> F {
        *index
    }

    fn chunk_bits(&self, M: usize) -> Vec<usize> {
        let log2_M = M.ilog2() as usize;
        let bound_bits = BOUND.ilog2() as usize;

        let remainder_bits = if bound_bits % log2_M != 0 {
            let reminder = 1 << (bound_bits % log2_M);
            let cutoff = reminder + BOUND % M as u64;
            let cutoff_log2 = cutoff.ilog2() as usize;
            vec![cutoff_log2]
        } else {
            vec![]
        };
        iter::repeat(log2_M)
            .take(bound_bits / log2_M)
            .chain(remainder_bits)
            .collect_vec()
    }

    fn subtable_indices(&self, index_bits: Vec<bool>, log_M: usize) -> Vec<Vec<bool>> {
        index_bits.chunks(log_M).map(Vec::from).collect_vec()
    }
}

#[derive(Clone, Debug, Default, Copy)]
pub struct SmallBoundLookup<const BOUND: u64>;

impl<const BOUND: u64> LookupType for SmallBoundLookup<BOUND> {
    fn combine_lookups<F: PrimeField>(&self, operands: &[F], _: usize, M: usize) -> F {
        let weight = F::from(M as u64);
        inner_product(
            operands,
            iter::successors(Some(F::ONE), |power_of_weight| {
                Some(*power_of_weight * weight)
            })
            .take(operands.len())
            .collect_vec()
            .iter(),
        )
    }

    fn combine_lookup_expressions<F: PrimeField, E: ExtensionField<F>>(
        &self,
        expressions: Vec<Expression<E, usize>>,
        C: usize,
        M: usize,
    ) -> Expression<E, usize> {
        if expressions.len() == 1 {
            return expressions[0].clone();
        }
        Expression::distribute_powers(expressions, E::from_bases(&[F::from(M as u64)]))
    }

    // SubtableIndices map subtable to memories
    fn subtables<F: PrimeField, E: ExtensionField<F>>(
        &self,
        C: usize,
        M: usize,
    ) -> Vec<(Box<dyn LassoSubtable<F, E>>, SubtableIndices)> {
        vec![(
            Box::new(BoundSubtable::<F, E, BOUND>::new()),
            SubtableIndices::from(0),
        )]
    }

    fn output<F: PrimeField>(&self, index: &F) -> F {
        *index
    }

    fn chunk_bits(&self, M: usize) -> Vec<usize> {
        let log2_M = M.ilog2() as usize;
        let bound_bits = BOUND.ilog2() as usize;

        let reminder = 1 << (bound_bits % log2_M);
        let cutoff = reminder + BOUND % M as u64;
        let cutoff_log2 = cutoff.ilog2() as usize;
        vec![cutoff_log2]
    }

    fn subtable_indices(&self, index_bits: Vec<bool>, log_M: usize) -> Vec<Vec<bool>> {
        index_bits.chunks(log_M).map(Vec::from).collect_vec()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use gkr::ff_ext::ff::Field;
    use gkr::util::dev::std_rng;
    use gkr::{poly::box_dense_poly, util::dev::seeded_std_rng};
    use goldilocks::Goldilocks;
    use rand::RngCore;

    use crate::lasso::LassoSubtable;

    use super::BoundSubtable;

    type F = Goldilocks;
    const LOG2_M: usize = 16;
    const M: usize = 1 << LOG2_M;

    #[test]
    fn full_subtable_mle_eval_correct() {
        let s = FullLimbSubtable::<F, F>::new();

        let poly = box_dense_poly::<F, F, _>(s.materialize(M));
        let num_vars = poly.num_vars();
        let r = (0..num_vars)
            .map(|_| F::random(seeded_std_rng()))
            .collect::<Vec<_>>();

        let full_poly_eval = poly.evaluate(&r);
        let func_eval = s.evaluate_mle(&r, M);
        let term_eval = s.evaluate_mle_expr(LOG2_M).evaluate(&r);

        assert_eq!(full_poly_eval, func_eval);
        assert_eq!(full_poly_eval, term_eval);
    }

    #[test]
    fn bound_subtable_mle_eval_correct() {
        const BOUND: u64 = (1 << 55) + 55;

        let s = BoundSubtable::<F, F, BOUND>::new();
        let evals = s.materialize(M);

        let poly = box_dense_poly::<F, F, _>(evals);
        let num_vars = poly.num_vars();
        let rng = &mut std_rng();
        let r = (0..num_vars)
            .map(|_| F::from(rng.next_u64()))
            .collect::<Vec<_>>();

        let full_poly_eval = poly.evaluate(&r);
        let func_eval = s.evaluate_mle(&r, M);
        let term_eval = s.evaluate_mle_expr(LOG2_M).evaluate(&r);

        assert_eq!(full_poly_eval, func_eval);
        assert_eq!(full_poly_eval, term_eval);
    }
}
