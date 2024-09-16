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

    fn subtable_id(&self) -> super::SubtableId {
        "full".to_string()
    }
}

impl<F, E> FullLimbSubtable<F, E> {
    pub fn new() -> Self {
        Self(PhantomData)
    }
}

#[derive(Clone, Debug, Default)]
pub struct BoundSubtable<F, E> {
    bound: u64,
    _marker: PhantomData<(F, E)>,
}

impl<F: PrimeField, E: ExtensionField<F>> LassoSubtable<F, E> for BoundSubtable<F, E> {
    fn materialize(&self, M: usize) -> Vec<F> {
        let bound_bits = self.bound.ilog2() as usize;
        let reminder = 1 << (bound_bits % M.ilog2() as usize);
        let cutoff = (reminder + self.bound % M as u64) as usize;

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

        let bound_bits = self.bound.ilog2() as usize;
        let reminder = 1 << (bound_bits % log2_M);
        let cutoff = reminder + self.bound % (1 << log2_M);
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
                        let mut term = E::from_bases(&[F::from(g_base + k)]);
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
        let bound_bits = self.bound.ilog2() as usize;
        let reminder = 1 << (bound_bits % log2_M);
        let cutoff = reminder + self.bound % M;
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
                    let mut term = PolyExpr::u64(g_base + k);
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

    fn subtable_id(&self) -> super::SubtableId {
        format!("bound_{}", self.bound)
    }
}

impl<F, E> BoundSubtable<F, E> {
    pub fn new(bound: u64) -> Self {
        Self {
            bound,
            _marker: PhantomData,
        }
    }
}

#[derive(Clone, Debug, Default, Copy)]
pub struct RangeLookup<F, E> {
    bound: u64,
    _marker: PhantomData<(F, E)>,
}

impl<F: PrimeField, E: ExtensionField<F>> LookupType<F, E> for RangeLookup<F, E> {
    fn combine_lookups(&self, operands: &[F], _: usize, M: usize) -> F {
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

    fn combine_lookup_expressions(
        &self,
        expressions: Vec<Expression<E, usize>>,
        _C: usize,
        M: usize,
    ) -> Expression<E, usize> {
        Expression::distribute_powers(expressions, E::from_bases(&[F::from(M as u64)]))
    }

    // SubtableIndices map subtable to memories
    fn subtables(
        &self,
        _C: usize,
        M: usize,
    ) -> Vec<(Box<dyn LassoSubtable<F, E>>, SubtableIndices)> {
        let full = Box::new(FullLimbSubtable::<F, E>::new());
        let log_M = M.ilog2() as usize;
        let bound_bits = self.bound.ilog2() as usize;
        let num_chunks = bound_bits / log_M;
        let rem = Box::new(BoundSubtable::<F, E>::new(self.bound));

        if self.bound % M as u64 == 0 {
            vec![(full, SubtableIndices::from(0..num_chunks))]
        } else if self.bound < M as u64 {
            vec![(rem, SubtableIndices::from(0))]
        } else {
            vec![
                (full, SubtableIndices::from(0..num_chunks)),
                (rem, SubtableIndices::from(num_chunks)),
            ]
        }
    }

    fn output(&self, index: &F) -> F {
        *index
    }

    fn chunk_bits(&self, M: usize) -> Vec<usize> {
        let log2_M = M.ilog2() as usize;
        let bound_bits = self.bound.ilog2() as usize;

        let remainder_bits = if self.bound % M as u64 != 0 {
            let reminder = 1 << (bound_bits % log2_M);
            let cutoff = reminder + self.bound % M as u64;
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

    fn lookup_id(&self) -> super::LookupId {
        format!("range_{}", self.bound)
    }
}

impl<F: PrimeField, E: ExtensionField<F>> RangeLookup<F, E> {
    pub fn new_boxed(bound: u64) -> Box<dyn LookupType<F, E>> {
        Box::new(Self {
            bound,
            _marker: PhantomData,
        })
    }
}

impl RangeLookup<(), ()> {
    pub fn id_for(bound: u64) -> super::LookupId {
        format!("range_{}", bound)
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

        let s = BoundSubtable::<F, F>::new(BOUND);
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
