use std::{iter, marker::PhantomData};

use gkr::ff_ext::{ff::PrimeField, ExtensionField};
use itertools::Itertools;

use gkr::{
    poly::{MultilinearPolyTerms, PolyExpr},
    util::{
        arithmetic::inner_product,
        expression::Expression,
    },
};

use gkr::circuit::node::lasso::{LassoSubtable, LookupType, SubtableIndices};

#[derive(Clone, Debug, Default)]
pub struct FullLimbSubtable<F, E>(PhantomData<(F, E)>);

impl<F: PrimeField, E: ExtensionField<F>> LassoSubtable<F, E> for FullLimbSubtable<F, E> {
    fn materialize(&self, M: usize) -> Vec<F> {
        (0..M).map(|x| F::from(x as u64)).collect_vec()
    }

    // fn evaluate_mle(&self, point: &[E]) -> E {
    //     let b = point.len();
    //     let mut result = E::ZERO;
    //     for i in 0..b {
    //         result += point[b] * F::from(1u64 << (i));
    //     }
    //     result
    // }

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
pub struct ReminderSubtable<F, E, const NUM_BITS: usize>(PhantomData<(F, E)>);

impl<F: PrimeField, E: ExtensionField<F>, const NUM_BITS: usize> LassoSubtable<F, E>
    for ReminderSubtable<F, E, NUM_BITS>
{
    fn materialize(&self, M: usize) -> Vec<F> {
        let cutoff = 1 << (NUM_BITS % M.ilog2() as usize);

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

    // fn evaluate_mle(&self, point: &[E]) -> E {
    //     let b = point.len();
    //     let remainder = todo!();
    //     // let remainder = NUM_BITS % LIMB_SIZE;
    //     let mut result = E::ZERO;
    //     for i in 0..b {
    //         if i < remainder {
    //             result += point[b] * F::from(1u64 << (i));
    //         } else {
    //             result *= E::ONE - point[b];
    //         }
    //     }
    //     result
    // }

    fn evaluate_mle_expr(&self, log2_M: usize) -> MultilinearPolyTerms<F> {
        let remainder = NUM_BITS % log2_M;
        let rem_init = PolyExpr::Var(0);
        let mut rem_terms = vec![rem_init];
        (1..remainder).for_each(|i| {
            let coeff = PolyExpr::Pow(Box::new(PolyExpr::Const(F::from(2))), i as u32);
            let x = PolyExpr::Var(i);
            let term = PolyExpr::Prod(vec![coeff, x]);
            rem_terms.push(term);
        });

        let mut other_terms = vec![];

        (remainder..log2_M).for_each(|i| {
            let x = PolyExpr::Var(i);
            let x_neg = PolyExpr::Prod(vec![PolyExpr::Const(F::ZERO - F::ONE), x]);
            let term = PolyExpr::Sum(vec![PolyExpr::Const(F::ONE), x_neg]);
            other_terms.push(term);
        });
        MultilinearPolyTerms::new(
            log2_M,
            PolyExpr::Prod([vec![PolyExpr::Sum(rem_terms)], other_terms].concat()),
        )
    }
}

impl<F, E, const NUM_BITS: usize> ReminderSubtable<F, E, NUM_BITS> {
    pub fn new() -> Self {
        Self(PhantomData)
    }
}

#[derive(Clone, Debug, Default, Copy)]
pub struct RangeStategy<const NUM_BITS: usize>;

impl<const NUM_BITS: usize> LookupType for RangeStategy<NUM_BITS> {
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
        let num_chunks = NUM_BITS / log_M;
        if NUM_BITS % log_M == 0 {
            vec![(full, SubtableIndices::from(0..num_chunks))]
        } else {
            let rem = Box::new(ReminderSubtable::<F, E, NUM_BITS>::new());
            vec![
                (full, SubtableIndices::from(0..num_chunks)),
                (rem, SubtableIndices::from(num_chunks)),
            ]
        }
    }

    fn output<F: PrimeField>(&self, index: &F) -> F {
        *index
    }

    fn chunk_bits(&self, log_M: usize) -> Vec<usize> {
        let remainder_bits = if NUM_BITS % log_M != 0 {
            vec![NUM_BITS % log_M]
        } else {
            vec![]
        };
        iter::repeat(log_M)
            .take(NUM_BITS / log_M)
            .chain(remainder_bits)
            .collect_vec()
    }

    // fn to_indices<F: PrimeField>(&self, value: &F) -> Vec<usize> {
    //     chunk_operand_usize(value.into(), div_ceil(NUM_BITS, LIMB_BITS), LIMB_BITS)
    // }

    fn subtable_indices(&self, index_bits: Vec<bool>, log_M: usize) -> Vec<Vec<bool>> {
        index_bits.chunks(log_M).map(Vec::from).collect_vec()
    }

    // fn num_memories(&self) -> usize {
    //     div_ceil(NUM_BITS, LIMB_BITS)
    // }
}

#[derive(Clone, Debug)]
pub struct RangeTable<F, E, const NUM_BITS: usize, const LIMB_BITS: usize>(
    PhantomData<F>,
    PhantomData<E>,
);

impl<F, E, const NUM_BITS: usize, const LIMB_BITS: usize> RangeTable<F, E, NUM_BITS, LIMB_BITS> {
    pub fn new() -> Self {
        Self(PhantomData, PhantomData)
    }
}

pub fn chunk_operand_usize(x: u64, C: usize, chunk_len: usize) -> Vec<usize> {
    let bit_mask = (1 << chunk_len) - 1;
    (0..C)
        .map(|i| {
            let shift = ((C - i - 1) * chunk_len) as u32;
            (x.checked_shr(shift).unwrap_or(0) & bit_mask) as usize
        })
        .collect()
}
