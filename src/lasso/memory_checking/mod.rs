pub mod prover;
pub mod verifier;

use gkr::{
    ff_ext::{ff::PrimeField, ExtensionField},
    poly::BoxMultilinearPoly,
};
use itertools::{chain, Itertools};
pub use prover::MemoryCheckingProver;

#[derive(Debug)]
struct MemoryGKR<'a, F: PrimeField, E: ExtensionField<F>> {
    init: BoxMultilinearPoly<'a, F, E>,
    read: BoxMultilinearPoly<'a, F, E>,
    write: BoxMultilinearPoly<'a, F, E>,
    final_read: BoxMultilinearPoly<'a, F, E>,
}

impl<'a, F: PrimeField, E: ExtensionField<F>> MemoryGKR<'a, F, E> {
    pub fn new(
        init: BoxMultilinearPoly<'a, F, E>,
        read: BoxMultilinearPoly<'a, F, E>,
        write: BoxMultilinearPoly<'a, F, E>,
        final_read: BoxMultilinearPoly<'a, F, E>,
    ) -> Self {
        Self {
            init,
            read,
            write,
            final_read,
        }
    }
}

type Poly<'a, F, E> = BoxMultilinearPoly<'a, F, E>;

#[derive(Clone, Debug)]
pub struct Chunk<'a, F: PrimeField, E: ExtensionField<F>> {
    pub(super) chunk_index: usize,
    pub(super) dim: &'a BoxMultilinearPoly<'a, F, E>,
    pub(super) read_ts_poly: &'a BoxMultilinearPoly<'a, F, E>,
    pub(super) final_cts_poly: &'a BoxMultilinearPoly<'a, F, E>,
    pub(super) memories: Vec<Memory<'a, F, E>>,
}

impl<'a, F: PrimeField, E: ExtensionField<F>> Chunk<'a, F, E> {
    pub(crate) fn new(
        chunk_index: usize,
        dim: &'a BoxMultilinearPoly<'a, F, E>,
        read_ts_poly: &'a BoxMultilinearPoly<'a, F, E>,
        final_cts_poly: &'a BoxMultilinearPoly<'a, F, E>,
        memory: Memory<'a, F, E>,
    ) -> Self {
        // sanity check
        assert_eq!(dim.num_vars(), read_ts_poly.num_vars());

        Self {
            chunk_index,
            dim,
            read_ts_poly,
            final_cts_poly,
            memories: vec![memory],
        }
    }

    pub fn chunk_index(&self) -> usize {
        self.chunk_index
    }

    pub fn chunk_bits(&self) -> usize {
        self.final_cts_poly.num_vars()
    }

    pub fn num_reads(&self) -> usize {
        1 << self.dim.num_vars()
    }

    pub fn chunk_polys(&self) -> impl Iterator<Item = &BoxMultilinearPoly<F, E>> {
        chain!([self.dim, self.read_ts_poly, self.final_cts_poly])
    }

    pub fn chunk_poly_evals(&self, x: &[E], y: &[E]) -> Vec<E> {
        vec![
            self.dim.evaluate(x),
            self.read_ts_poly.evaluate(x),
            self.final_cts_poly.evaluate(y),
        ]
    }

    pub fn e_poly_evals(&self, x: &[E]) -> Vec<E> {
        self.memories
            .iter()
            .map(|memory| memory.e_poly.evaluate(x))
            .collect_vec()
    }

    pub(super) fn memories(&self) -> impl Iterator<Item = &Memory<F, E>> {
        self.memories.iter()
    }

    pub(super) fn add_memory(&mut self, memory: Memory<'a, F, E>) {
        // sanity check
        let chunk_bits = self.chunk_bits();
        let num_reads = self.num_reads();
        assert_eq!(chunk_bits, memory.subtable_poly.num_vars());
        println!(
            "num_reads: {}, memory.e_poly.num_vars(): {}",
            num_reads,
            memory.e_poly.num_vars()
        );
        assert_eq!(num_reads, 1 << memory.e_poly.num_vars());

        self.memories.push(memory);
    }
}

#[derive(Clone, Debug)]
pub(super) struct Memory<'a, F: PrimeField, E: ExtensionField<F>> {
    subtable_poly: &'a BoxMultilinearPoly<'a, F, E>,
    pub(crate) e_poly: &'a BoxMultilinearPoly<'a, F, E>,
}

impl<'a, F: PrimeField, E: ExtensionField<F>> Memory<'a, F, E> {
    pub(crate) fn new(
        subtable_poly: &'a BoxMultilinearPoly<'a, F, E>,
        e_poly: &'a BoxMultilinearPoly<'a, F, E>,
    ) -> Self {
        Self {
            subtable_poly,
            e_poly,
        }
    }

    pub fn polys(&'a self) -> impl Iterator<Item = &'a BoxMultilinearPoly<'a, F, E>> {
        chain!([self.subtable_poly, self.e_poly])
    }
}
