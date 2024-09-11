use fixedbitset::FixedBitSet;
use std::ops::Range;
use std::{any::TypeId, fmt::Debug};
use strum::{EnumCount, IntoEnumIterator};
use enum_dispatch::enum_dispatch;
use gkr::{
    ff_ext::{ff::PrimeField, ExtensionField},
    poly::MultilinearPolyTerms,
    util::expression::Expression,
};

pub mod range;
use range::*;
use crate::sk_encryption_circuit::{RangeLookups, RangeSubtables, BfvEncryptBlock};

pub type SubtableId = TypeId;
pub type LookupId = TypeId;

#[enum_dispatch]
pub trait LassoSubtable<F: PrimeField, E: ExtensionField<F>>: 'static + Sync + Debug {
    /// Returns the TypeId of this subtable.
    /// The `Jolt` trait has associated enum types `InstructionSet` and `Subtables`.
    /// This function is used to resolve the many-to-many mapping between `InstructionSet` variants
    /// and `Subtables` variants,
    fn subtable_id(&self) -> SubtableId {
        TypeId::of::<Self>()
    }
    /// Fully materializes a subtable of size `M`, reprensented as a Vec of length `M`.
    fn materialize(&self, M: usize) -> Vec<F>;

    // fn evaluate_mle(&self, point: &[E]) -> E;

    /// Expression to evaluate the multilinear extension polynomial for this subtable at the given `point`,
    /// interpreted to be of size log_2(M), where M is the size of the subtable.
    fn evaluate_mle_expr(&self, log2_M: usize) -> MultilinearPolyTerms<F>;
}

pub trait SubtableSet<F: PrimeField, E: ExtensionField<F>>:
    LassoSubtable<F, E> + IntoEnumIterator + EnumCount + From<SubtableId> + Send + Sync + Debug
{
    fn enum_index(subtable: Box<dyn LassoSubtable<F, E>>) -> usize {
        let s = Self::from(subtable.subtable_id());
        let byte = unsafe { *(&s as *const Self as *const u8) };
        byte as usize
    }
}

pub trait CircuitLookups:
    LookupType + IntoEnumIterator + EnumCount + Send + Sync + std::fmt::Debug + Copy
{
    fn enum_index(lookup_type: &Self) -> usize {
        // Discriminant: https://doc.rust-lang.org/reference/items/enumerations.html#pointer-casting
        let byte = unsafe { *(lookup_type as *const Self as *const u8) };
        byte as usize
    }
}

#[enum_dispatch]
pub trait LookupType: Clone + Send + Sync {
    // /// Returns the TypeId of this lookup type.
    // fn lookup_id(&self) -> LookupId {
    //     TypeId::of::<Self>()
    // }

    /// The `g` function that computes T[r] = g(T_1[r_1], ..., T_k[r_1], T_{k+1}[r_2], ..., T_{\alpha}[r_c])
    fn combine_lookups<F: PrimeField>(&self, operands: &[F], C: usize, M: usize) -> F;

    fn combine_lookup_expressions<F: PrimeField, E: ExtensionField<F>>(
        &self,
        expressions: Vec<Expression<E, usize>>,
        C: usize,
        M: usize,
    ) -> Expression<E, usize>;

    /// Returns a Vec of the unique subtable types used by this instruction. For some instructions,
    /// e.g. SLL, the list of subtables depends on the dimension `C`.
    fn subtables<F: PrimeField, E: ExtensionField<F>>(
        &self,
        C: usize,
        M: usize,
    ) -> Vec<(Box<dyn LassoSubtable<F, E>>, SubtableIndices)>;

    // fn to_indices<F: PrimeField>(&self, value: &F) -> Vec<usize>;

    fn output<F: PrimeField>(&self, index: &F) -> F;

    fn chunk_bits(&self, log_M: usize) -> Vec<usize>;

    /// Returns the indices of each subtable lookups
    /// The length of `index_bits` is same as actual bit length of table index
    fn subtable_indices(&self, index_bits: Vec<bool>, log_M: usize) -> Vec<Vec<bool>>;

    // fn num_memories(&self) -> usize;
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

#[macro_export]
macro_rules! subtable_enum {
    ($enum_name:ident, $($alias:ident: $struct:ty),+) => {
        #[enum_dispatch(LassoSubtable<F, E>)]
        #[derive(EnumCount, EnumIter, Debug)]
        pub enum $enum_name<F: PrimeField, E: ExtensionField<F>> { $($alias($struct)),+ }
        impl<F: PrimeField, E: ExtensionField<F>> From<SubtableId> for $enum_name<F, E> {
          fn from(subtable_id: SubtableId) -> Self {
            $(
              if subtable_id == TypeId::of::<$struct>() {
                $enum_name::from(<$struct>::new())
              } else
            )+
            { panic!("Unexpected subtable id {:?}", subtable_id) }
          }
        }

        impl<F: PrimeField, E: ExtensionField<F>> SubtableSet<F, E> for $enum_name<F, E> {}
    };
}
