use gkr::{
    ff_ext::{
        ff::{Field, PrimeField},
        ExtensionField,
    },
    sum_check::verify_sum_check,
    transcript::TranscriptRead,
    util::arithmetic::inner_product,
    Error,
};
use itertools::{izip, Itertools};
use std::{iter, marker::PhantomData};

use crate::lasso::LassoSubtable;

use super::MemoryCheckingProver;

#[derive(Debug)]
pub struct Chunk<F: Field, E: Field, const M: usize> {
    chunk_index: usize,
    chunk_bits: usize,
    pub(crate) memory: Vec<Memory<F, E>>,
}

impl<F: PrimeField, E: ExtensionField<F>, const M: usize> Chunk<F, E, M> {
    pub fn chunk_polys_index(&self, offset: usize, num_chunks: usize) -> Vec<usize> {
        let dim_poly_index = offset + 1 + self.chunk_index;
        let read_ts_poly_index = offset + 1 + num_chunks + self.chunk_index;
        let final_cts_poly_index = offset + 1 + 2 * num_chunks + self.chunk_index;
        vec![dim_poly_index, read_ts_poly_index, final_cts_poly_index]
    }

    pub fn new(chunk_index: usize, chunk_bits: usize, memory: Memory<F, E>) -> Self {
        Self {
            chunk_index,
            chunk_bits,
            memory: vec![memory],
        }
    }

    pub fn num_memories(&self) -> usize {
        self.memory.len()
    }

    pub fn chunk_bits(&self) -> usize {
        self.chunk_bits
    }

    pub fn add_memory(&mut self, memory: Memory<F, E>) {
        self.memory.push(memory);
    }

    pub fn memory_indices(&self) -> Vec<usize> {
        self.memory
            .iter()
            .map(|memory| memory.memory_index)
            .collect_vec()
    }

    #[allow(clippy::too_many_arguments)]
    pub fn verify_memories(
        &self,
        read_xs: &[E],
        write_xs: &[E],
        init_ys: &[E],
        final_read_ys: &[E],
        y: &[E],
        hash: impl Fn(&E, &E, &E) -> E,
        transcript: &mut dyn TranscriptRead<F, E>,
    ) -> Result<(), Error> {
        let [dim_x, read_ts_poly_x, final_cts_poly_y] =
            transcript.read_felt_exts(3)?.try_into().unwrap();
        let e_poly_xs = transcript.read_felt_exts(self.num_memories())?;
        let id_poly_y = inner_product(
            iter::successors(Some(E::ONE), |power_of_two| Some(power_of_two.double()))
                .take(y.len()),
            y.to_vec(),
        );
        self.memory.iter().enumerate().for_each(|(i, memory)| {
            assert_eq!(read_xs[i], hash(&dim_x, &e_poly_xs[i], &read_ts_poly_x));
            assert_eq!(
                write_xs[i],
                hash(&dim_x, &e_poly_xs[i], &(read_ts_poly_x + E::ONE))
            );
            let subtable_poly_y = memory.subtable.evaluate_mle(y, M);

            assert_eq!(init_ys[i], hash(&id_poly_y, &subtable_poly_y, &E::ZERO));
            assert_eq!(
                final_read_ys[i],
                hash(&id_poly_y, &subtable_poly_y, &final_cts_poly_y)
            );
        });
        // Ok((dim_x, read_ts_poly_x, final_cts_poly_y, e_poly_xs))
        Ok(())
    }
}

#[derive(Debug)]
pub struct Memory<F: Field, E: Field> {
    memory_index: usize,
    subtable: Box<dyn LassoSubtable<F, E>>,
}

impl<F: Field, E: Field> Memory<F, E> {
    pub fn new(memory_index: usize, subtable: Box<dyn LassoSubtable<F, E>>) -> Self {
        Self {
            memory_index,
            subtable,
        }
    }
}

#[derive(Debug)]
pub struct MemoryCheckingVerifier<F: PrimeField, E: ExtensionField<F>, const M: usize> {
    /// chunks with the same bits size
    chunks: Vec<Chunk<F, E, M>>,
    _marker: PhantomData<F>,
    _marker_e: PhantomData<E>,
}

impl<F: PrimeField, E: ExtensionField<F>, const M: usize> MemoryCheckingVerifier<F, E, M> {
    pub fn new(chunks: Vec<Chunk<F, E, M>>) -> Self {
        Self {
            chunks,
            _marker: PhantomData,
            _marker_e: PhantomData,
        }
    }

    pub fn verify(
        &self,
        num_reads: usize,
        gamma: &E,
        tau: &E,
        transcript: &mut dyn TranscriptRead<F, E>,
    ) -> Result<(), Error> {
        // warning: this cast is insecure
        // TODO: need to rewrite memory checking to work over the extension field
        let tau = tau.as_bases()[0];
        let gamma = gamma.as_bases()[0];

        let num_memories: usize = self.chunks.iter().map(|chunk| chunk.num_memories()).sum();
        let memory_bits = self.chunks[0].chunk_bits();
        let (read_write_xs, _) = Self::verify_grand_product(
            num_reads,
            iter::repeat(None).take(2 * num_memories),
            transcript,
        )?;
        let (read_xs, write_xs) = read_write_xs.split_at(num_memories);

        let (init_final_read_ys, y) = Self::verify_grand_product(
            memory_bits,
            iter::repeat(None).take(2 * num_memories),
            transcript,
        )?;
        let (init_ys, final_read_ys) = init_final_read_ys.split_at(num_memories);

        let hash = |a: &E, v: &E, t: &E| -> E { *a + *v * gamma + *t * gamma.square() - tau };
        let mut offset = 0;

        self.chunks.iter().for_each(|chunk| {
            let num_memories = chunk.num_memories();
            let _ = chunk.verify_memories(
                &read_xs[offset..offset + num_memories],
                &write_xs[offset..offset + num_memories],
                &init_ys[offset..offset + num_memories],
                &final_read_ys[offset..offset + num_memories],
                &y,
                hash,
                transcript,
            );
            offset += num_memories;
        });

        Ok(())
    }

    fn verify_grand_product(
        num_vars: usize,
        claimed_v_0s: impl IntoIterator<Item = Option<E>>,
        transcript: &mut dyn TranscriptRead<F, E>,
    ) -> Result<(Vec<E>, Vec<E>), Error> {
        let claimed_v_0s = claimed_v_0s.into_iter().collect_vec();
        let num_batching = claimed_v_0s.len();

        assert!(num_batching != 0);
        let claimed_v_0s = {
            claimed_v_0s
                .into_iter()
                .map(|claimed| match claimed {
                    Some(claimed) => {
                        transcript.common_felts(claimed.as_bases());
                        Ok(claimed)
                    }
                    None => transcript.read_felt_ext(),
                })
                .try_collect::<_, Vec<_>, _>()?
        };

        (0..num_vars).try_fold((claimed_v_0s, Vec::new()), |result, num_vars| {
            let (claimed_v_ys, _) = result;

            let (mut x, evals) = if num_vars == 0 {
                let evals = transcript.read_felt_exts(2 * num_batching)?;
                for (claimed_v, (&v_l, &v_r)) in izip!(claimed_v_ys, evals.iter().tuples()) {
                    if claimed_v != v_l * v_r {
                        return Err(Error::InvalidSumCheck(
                            "unmatched sum check output".to_string(),
                        ));
                    }
                }

                (Vec::new(), evals)
            } else {
                let gamma = transcript.squeeze_challenge();
                let g = MemoryCheckingProver::sum_check_function(num_vars, num_batching, gamma);

                let (_x_eval, x) = {
                    let claim = MemoryCheckingProver::sum_check_claim(&claimed_v_ys, gamma);
                    verify_sum_check(&g, claim, transcript)?
                };

                let evals = transcript.read_felt_exts(2 * num_batching)?;

                (x, evals)
            };

            let mu = transcript.squeeze_challenge();

            let v_xs = MemoryCheckingProver::layer_down_claim(&evals, mu);
            x.push(mu);

            Ok((v_xs, x))
        })
    }
}
