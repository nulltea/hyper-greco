use enum_dispatch::enum_dispatch;
use fixedbitset::FixedBitSet;
use gkr::{
    ff_ext::{ff::PrimeField, ExtensionField},
    poly::MultilinearPolyTerms,
    util::expression::Expression,
};
use std::ops::Range;
use std::{any::TypeId, fmt::Debug};
use strum::{EnumCount, IntoEnumIterator};

pub mod range;
use range::*;
// for some reason #[enum_dispatch] needs this
use crate::sk_encryption_circuit::*;

pub type SubtableId = String;
pub type LookupId = String;

pub trait LassoSubtable<F: PrimeField, E: ExtensionField<F>>: 'static + Sync + Debug {
    /// Returns the TypeId of this subtable.
    /// The `Jolt` trait has associated enum types `InstructionSet` and `Subtables`.
    /// This function is used to resolve the many-to-many mapping between `InstructionSet` variants
    /// and `Subtables` variants,
    fn subtable_id(&self) -> SubtableId;

    /// Fully materializes a subtable of size `M`, reprensented as a Vec of length `M`.
    fn materialize(&self, M: usize) -> Vec<F>;

    fn evaluate_mle(&self, point: &[E], M: usize) -> E;

    /// Expression to evaluate the multilinear extension polynomial for this subtable at the given `point`,
    /// interpreted to be of size log_2(M), where M is the size of the subtable.
    fn evaluate_mle_expr(&self, log2_M: usize) -> MultilinearPolyTerms<F>;
}

pub trait LookupType<F: PrimeField, E: ExtensionField<F>>: 'static + Send + Sync + Debug + LookupClone<F, E> {
    /// Returns the identifier of this lookup type.
    fn lookup_id(&self) -> LookupId;

    /// The `g` function that computes T[r] = g(T_1[r_1], ..., T_k[r_1], T_{k+1}[r_2], ..., T_{\alpha}[r_c])
    fn combine_lookups(&self, operands: &[F], C: usize, M: usize) -> F;

    fn combine_lookup_expressions(
        &self,
        expressions: Vec<Expression<E, usize>>,
        C: usize,
        M: usize,
    ) -> Expression<E, usize>;

    /// Returns a Vec of the unique subtable types used by this instruction. For some instructions,
    /// e.g. SLL, the list of subtables depends on the dimension `C`.
    fn subtables(&self, C: usize, M: usize)
        -> Vec<(Box<dyn LassoSubtable<F, E>>, SubtableIndices)>;

    // fn to_indices<F: PrimeField>(&self, value: &F) -> Vec<usize>;

    fn output(&self, index: &F) -> F;

    fn chunk_bits(&self, M: usize) -> Vec<usize>;

    /// Returns the indices of each subtable lookups
    /// The length of `index_bits` is same as actual bit length of table index
    fn subtable_indices(&self, index_bits: Vec<bool>, log_M: usize) -> Vec<Vec<bool>>;

    // fn num_memories(&self) -> usize;
}

pub trait LookupClone<F, E> {
    fn clone_box(&self) -> Box<dyn LookupType<F, E>>;
}

impl<T, F: PrimeField, E: ExtensionField<F>> LookupClone<F, E> for T
where
    T: LookupType<F, E> + Clone + 'static,
{
    fn clone_box(&self) -> Box<dyn LookupType<F, E>> {
        Box::new(self.clone())
    }
}

impl<F, E> Clone for Box<dyn LookupType<F, E>> {
    fn clone(&self) -> Self {
        self.clone_box()
    }
}

#[derive(Clone)]
pub struct SubtableIndices {
    bitset: FixedBitSet,
}

impl SubtableIndices {
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            bitset: FixedBitSet::with_capacity(capacity),
        }
    }

    pub fn union_with(&mut self, other: &Self) {
        self.bitset.union_with(&other.bitset);
    }

    pub fn iter(&self) -> impl Iterator<Item = usize> + '_ {
        self.bitset.ones()
    }

    #[allow(clippy::len_without_is_empty)]
    pub fn len(&self) -> usize {
        self.bitset.count_ones(..)
    }

    pub fn contains(&self, index: usize) -> bool {
        self.bitset.contains(index)
    }
}

impl From<usize> for SubtableIndices {
    fn from(index: usize) -> Self {
        let mut bitset = FixedBitSet::new();
        bitset.grow_and_insert(index);
        Self { bitset }
    }
}

impl From<Range<usize>> for SubtableIndices {
    fn from(range: Range<usize>) -> Self {
        let bitset = FixedBitSet::from_iter(range);
        Self { bitset }
    }
}
