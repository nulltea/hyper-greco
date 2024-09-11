use itertools::{chain, izip, Itertools};
use plonkish_backend::util::parallel::{num_threads, parallelize_iter};
use rayon::prelude::{IntoParallelIterator, ParallelIterator};
use std::{array, iter};

use gkr::{
    ff_ext::{
        ff::{Field, PrimeField},
        ExtensionField,
    },
    poly::{box_dense_poly, BoxMultilinearPoly},
    sum_check::{generic::Generic, prove_sum_check, SumCheckPoly},
    transcript::TranscriptWrite,
    util::{
        arithmetic::{div_ceil, inner_product, powers},
        expression::Expression,
    },
    Error,
};

use super::{Chunk, MemoryGKR};

pub struct MemoryCheckingProver<'a, F: PrimeField, E: ExtensionField<F>> {
    /// chunks with the same bits size
    chunks: Vec<Chunk<'a, F, E>>,
    /// GKR initial polynomials for each memory
    pub memories: Vec<MemoryGKR<'a, F, E>>,
    // gamma: F,
}

impl<'a, F: PrimeField + Field, E: ExtensionField<F>> MemoryCheckingProver<'a, F, E> {
    // T_1[dim_1(x)], ..., T_k[dim_1(x)],
    // ...
    // T_{\alpha-k+1}[dim_c(x)], ..., T_{\alpha}[dim_c(x)]
    pub fn new(chunks: Vec<Chunk<'a, F, E>>, tau: &E, gamma: &E) -> Self {
        // TODO: suppoer extension degree > 1
        let tau = tau.as_bases()[0];
        let gamma = gamma.as_bases()[0];
        let num_reads = chunks[0].num_reads();
        println!("num_reads: {}", num_reads);
        let memory_size = 1 << chunks[0].chunk_bits();
        println!("memory_size: {}", memory_size);

        let hash = |a: &F, v: &F, t: &F| -> F { *a + *v * gamma + *t * gamma.square() - tau };

        let memories_gkr: Vec<MemoryGKR<F, E>> = (0..chunks.len())
            .into_par_iter()
            .flat_map(|i| {
                let chunk = &chunks[i];
                let chunk_polys = chunk.chunk_polys().collect_vec();
                let (dim, read_ts_poly, final_cts_poly) =
                    (chunk_polys[0], chunk_polys[1], chunk_polys[2]);
                chunk
                    .memories()
                    .map(|memory| {
                        let memory_polys = memory.polys().collect_vec();
                        let (subtable_poly, e_poly) = (memory_polys[0], memory_polys[1]);
                        let mut init = vec![];
                        let mut read = vec![];
                        let mut write = vec![];
                        let mut final_read = vec![];
                        (0..memory_size).for_each(|i| {
                            init.push(hash(&F::from(i as u64), &subtable_poly[i], &F::ZERO));
                            final_read.push(hash(
                                &F::from(i as u64),
                                &subtable_poly[i],
                                &final_cts_poly[i],
                            ));
                        });
                        (0..num_reads).for_each(|i| {
                            read.push(hash(&dim[i], &e_poly[i], &read_ts_poly[i]));
                            write.push(hash(&dim[i], &e_poly[i], &(read_ts_poly[i] + F::ONE)));
                        });
                        MemoryGKR::new(
                            box_dense_poly(init),
                            box_dense_poly(read),
                            box_dense_poly(write),
                            box_dense_poly(final_read),
                        )
                    })
                    .collect_vec()
            })
            .collect();

        Self {
            chunks,
            memories: memories_gkr,
            // gamma: *gamma,
        }
    }

    pub fn inits(&self) -> impl Iterator<Item = &BoxMultilinearPoly<'a, F, E>> {
        self.memories.iter().map(|memory| &memory.init)
    }

    pub fn reads(&self) -> impl Iterator<Item = &BoxMultilinearPoly<'a, F, E>> {
        self.memories.iter().map(|memory| &memory.read)
    }

    pub fn writes(&self) -> impl Iterator<Item = &BoxMultilinearPoly<'a, F, E>> {
        self.memories.iter().map(|memory| &memory.write)
    }

    pub fn final_reads(&self) -> impl Iterator<Item = &BoxMultilinearPoly<'a, F, E>> {
        self.memories.iter().map(|memory| &memory.final_read)
    }

    pub fn iter(
        &self,
    ) -> impl Iterator<
        Item = (
            &BoxMultilinearPoly<F, E>,
            &BoxMultilinearPoly<F, E>,
            &BoxMultilinearPoly<F, E>,
            &BoxMultilinearPoly<F, E>,
        ),
    > {
        self.memories.iter().map(|memory| {
            (
                &memory.init,
                &memory.read,
                &memory.write,
                &memory.final_read,
            )
        })
    }

    pub fn claimed_v_0s(&self) -> impl IntoIterator<Item = Vec<Option<F>>> {
        let (claimed_read_0s, claimed_write_0s, claimed_init_0s, claimed_final_read_0s) = self
            .iter()
            .map(|(init, read, write, final_read)| {
                let claimed_init_0 = init.to_dense().iter().product();
                let claimed_read_0 = read.to_dense().iter().product();
                let claimed_write_0 = write.to_dense().iter().product();
                let claimed_final_read_0 = final_read.to_dense().iter().product();

                // sanity check
                debug_assert_eq!(
                    claimed_init_0 * claimed_write_0,
                    claimed_read_0 * claimed_final_read_0,
                    "Multiset hashes don't match",
                );
                (
                    Some(claimed_read_0),
                    Some(claimed_write_0),
                    Some(claimed_init_0),
                    Some(claimed_final_read_0),
                )
            })
            .multiunzip::<(Vec<_>, Vec<_>, Vec<_>, Vec<_>)>();
        chain!([
            claimed_read_0s,
            claimed_write_0s,
            claimed_init_0s,
            claimed_final_read_0s
        ])
    }

    pub fn prove(
        &mut self,
        // points_offset: usize,
        // lookup_opening_points: &mut Vec<Vec<F>>,
        // lookup_opening_evals: &mut Vec<Evaluation<F>>,
        transcript: &mut dyn TranscriptWrite<F, E>,
    ) -> Result<(), Error> {
        let num_batching = self.memories.len() * 2;

        let (_, x) = Self::prove_grand_product(
            num_batching,
            chain!(self.reads(), self.writes()),
            transcript,
        )?;

        let (_, y) = Self::prove_grand_product(
            num_batching,
            chain!(self.inits(), self.final_reads()),
            transcript,
        )?;

        self.chunks.iter().for_each(|chunk| {
            let chunk_poly_evals = chunk.chunk_poly_evals(&x, &y);
            let e_poly_xs = chunk.e_poly_evals(&x);
            transcript.write_felt_exts(&chunk_poly_evals).unwrap();
            transcript.write_felt_exts(&e_poly_xs).unwrap();
        });

        Ok(())
    }

    fn prove_grand_product<'b>(
        num_batching: usize,
        vs: impl Iterator<Item = &'b BoxMultilinearPoly<'a, F, E>>,
        transcript: &mut dyn TranscriptWrite<F, E>,
    ) -> Result<(Vec<E>, Vec<E>), Error>
    where
        'a: 'b,
    {
        let bottom_layers = vs.map(Layer::bottom).collect_vec();
        let layers = iter::successors(bottom_layers.into(), |layers| {
            (layers[0].num_vars() > 0).then(|| layers.iter().map(Layer::up).collect())
        })
        .collect_vec();

        let claimed_v_0s = {
            let v_0s = chain![layers.last().unwrap()]
                .map(|layer| {
                    let [v_l, v_r] = layer.polys().map(|poly| poly[0]);
                    E::from_bases(&[v_l * v_r])
                })
                .collect_vec();

            let mut hash_to_transcript = |claimed: Vec<Option<E>>, computed: Vec<_>| {
                izip!(claimed, computed)
                    .map(|(claimed, computed)| match claimed {
                        Some(claimed) => {
                            if cfg!(feature = "sanity-check") {
                                assert_eq!(claimed, computed)
                            }
                            transcript.common_felts(&computed.as_bases());
                            Ok(computed)
                        }
                        None => transcript.write_felt_ext(&computed).map(|_| computed),
                    })
                    .try_collect::<_, Vec<_>, _>()
            };

            hash_to_transcript(iter::repeat(None).take(num_batching).collect_vec(), v_0s)?
        };

        layers
            .iter()
            .rev()
            .try_fold((claimed_v_0s, Vec::new()), |result, layers| {
                let (claimed_v_ys, _y) = result;

                let num_vars = layers[0].num_vars();
                let polys = layers.iter().flat_map(|layer| layer.polys());

                let (mut x, evals) = if num_vars == 0 {
                    (
                        vec![],
                        polys.map(|poly| E::from_bases(&[poly[0]])).collect_vec(),
                    )
                } else {
                    let gamma = transcript.squeeze_challenge();

                    let g = Self::sum_check_function(num_vars, num_batching, gamma);

                    let (_, x, evals) = {
                        let claim = Self::sum_check_claim(&claimed_v_ys, gamma);
                        let polys = polys
                            .map(|e_poly| {
                                SumCheckPoly::Base::<_, _, _, BoxMultilinearPoly<E, E>>(e_poly)
                            })
                            .collect_vec();

                        prove_sum_check(&g, claim, polys, transcript)?
                    };

                    (x, evals)
                };

                transcript.write_felt_exts(&evals)?;

                let mu = transcript.squeeze_challenge();

                let v_xs = Self::layer_down_claim(&evals, mu);
                x.push(mu);

                Ok((v_xs, x))
            })
    }

    pub fn sum_check_function(num_vars: usize, num_batching: usize, gamma: E) -> Generic<F, E> {
        let exprs = &(0..2 * num_batching)
            .map(|idx| Expression::poly(idx))
            .tuples()
            .map(|(ref v_l, ref v_r)| v_l * v_r)
            .collect_vec();
        let expr = Expression::distribute_powers(exprs, gamma);

        let eq_r_x = Expression::poly(0);

        Generic::new(num_vars, &(eq_r_x * expr))
    }

    pub fn sum_check_claim(claimed_v_ys: &[E], gamma: E) -> E {
        inner_product(
            claimed_v_ys.to_vec(),
            &powers(gamma).take(claimed_v_ys.len()).collect_vec(),
        )
    }

    pub fn layer_down_claim(evals: &[E], mu: E) -> Vec<E> {
        evals
            .iter()
            .tuples()
            .map(|(&v_l, &v_r)| v_l + mu * (v_r - v_l))
            .collect_vec()
    }
}

struct Layer<'a, F, E> {
    v_l: BoxMultilinearPoly<'a, F, E>,
    v_r: BoxMultilinearPoly<'a, F, E>,
}

impl<'a, F: PrimeField, E: ExtensionField<F>> From<[Vec<F>; 2]> for Layer<'a, F, E> {
    fn from(values: [Vec<F>; 2]) -> Self {
        let [v_l, v_r] = values.map(box_dense_poly);
        Self { v_l, v_r }
    }
}

impl<'a, F: PrimeField, E: ExtensionField<F>> Layer<'a, F, E> {
    fn bottom<'b>(v: &BoxMultilinearPoly<'a, F, E>) -> Layer<'b, F, E> {
        let mid = v.to_dense().len() >> 1;
        [&v.to_dense()[..mid], &v.to_dense()[mid..]]
            .map(ToOwned::to_owned)
            .into()
    }

    fn num_vars(&self) -> usize {
        self.v_l.num_vars()
    }

    fn polys(&self) -> [&BoxMultilinearPoly<'a, F, E>; 2] {
        [&self.v_l, &self.v_r]
    }

    fn poly_chunks(&self, chunk_size: usize) -> impl Iterator<Item = (&[F], &[F])> {
        let [v_l, v_r] = self
            .polys()
            .map(|poly| poly.as_dense().unwrap().chunks(chunk_size));
        izip!(v_l, v_r)
    }

    fn up<'b>(&self) -> Layer<'b, F, E> {
        assert!(self.num_vars() != 0);

        let len = 1 << self.num_vars();
        let chunk_size = div_ceil(len, num_threads()).next_power_of_two();

        let mut outputs: [_; 2] = array::from_fn(|_| vec![F::ZERO; len >> 1]);
        let (v_up_l, v_up_r) = outputs.split_at_mut(1);

        parallelize_iter(
            izip!(
                chain![v_up_l, v_up_r].flat_map(|v_up| v_up.chunks_mut(chunk_size)),
                self.poly_chunks(chunk_size),
            ),
            |(v_up, (v_l, v_r))| {
                izip!(v_up, v_l, v_r).for_each(|(v_up, v_l, v_r)| {
                    *v_up = *v_l * *v_r;
                })
            },
        );

        outputs.into()
    }
}
