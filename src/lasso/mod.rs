use ark_std::log2;
use gkr::{
    circuit::node::{CombinedEvalClaim, EvalClaim, Node},
    ff_ext::ff::PrimeField,
    izip_par,
    poly::{box_dense_poly, BoxMultilinearPoly, DensePolynomial, MultilinearPoly},
    sum_check::{generic::Generic, prove_sum_check, verify_sum_check, SumCheckPoly},
    transcript::{TranscriptRead, TranscriptWrite},
    util::{
        arithmetic::{ExtensionField, Field},
        chain,
        expression::Expression,
        izip, Itertools,
    },
    Error,
};
use memory_checking::{Chunk, Memory, MemoryCheckingProver};
use plonkish_backend::{
    poly::multilinear::MultilinearPolynomial,
    util::arithmetic::{fe_to_bits_le, usize_from_bits_le},
};
use rayon::prelude::*;
use std::{
    collections::{BTreeMap, HashMap},
    iter,
    sync::Arc,
};
use table::LookupId;
pub use table::{LassoSubtable, LookupType, SubtableId, SubtableIndices};
use tracing::info_span;

pub mod memory_checking;
pub mod table;

#[derive(Debug)]
pub struct LassoNode<F: Field, E, const C: usize, const M: usize> {
    num_vars: usize,
    preprocessing: LassoPreprocessing<F, E>,
    // inputs_arity: usize,
    lookups: Vec<LookupId>,
}

impl<F: PrimeField, E: ExtensionField<F>, const C: usize, const M: usize> Node<F, E>
    for LassoNode<F, E, C, M>
{
    fn is_input(&self) -> bool {
        false
    }

    fn log2_input_size(&self) -> usize {
        self.num_vars.max(log2(M) as usize)
    }

    fn log2_output_size(&self) -> usize {
        0
    }

    fn evaluate(&self, _: Vec<&BoxMultilinearPoly<F, E>>) -> BoxMultilinearPoly<'static, F, E> {
        box_dense_poly([F::ZERO])
    }

    #[tracing::instrument(skip_all, name = "LassoNode::prove_claim_reduction")]
    fn prove_claim_reduction(
        &self,
        _: CombinedEvalClaim<E>,
        inputs: Vec<&BoxMultilinearPoly<F, E>>,
        transcript: &mut dyn TranscriptWrite<F, E>,
    ) -> Result<Vec<Vec<EvalClaim<E>>>, Error> {
        // println!("combined claim {} {:?}", claim.points.len(), claim.value);

        // VanillaNode::new(input_arity, log2_sub_input_size, gates, num_reps)

        // let [input, output] = inputs.try_into().unwrap();
        let polys = self.polynomialize(inputs[0]);
        let mock_lookup = self.preprocessing.lookups.values().next().unwrap();

        let LassoPolynomials {
            dims,
            read_cts: read_ts_polys,
            final_cts: final_cts_polys,
            e_polys,
            lookup_outputs,
            lookup_flag_polys,
            lookup_flag_bitvectors: _,
        } = polys;

        assert!(inputs[0].to_dense() == lookup_outputs.to_dense());

        // let [e_polys, dims, read_ts_polys, final_cts_polys] = polys;

        let num_vars = lookup_outputs.num_vars();
        assert_eq!(num_vars, self.num_vars);

        assert_eq!(final_cts_polys[0].num_vars(), log2(M) as usize);

        // should this be passed from CombinedEvalClaim?
        let r = transcript.squeeze_challenges(num_vars);

        let res = self.prove_collation_sum_check(
            &lookup_outputs,
            mock_lookup,
            &e_polys,
            &lookup_flag_polys,
            &r,
            num_vars,
            transcript,
        )?;

        let lookup_output_eval_claim = res.into_iter().take(1).collect_vec();

        let [gamma, tau] = transcript.squeeze_challenges(2).try_into().unwrap();

        // memory_checking
        Self::prove_memory_checking(
            &self.preprocessing,
            &dims,
            &read_ts_polys,
            &final_cts_polys,
            &e_polys,
            &gamma,
            &tau,
            transcript,
        )?;

        Ok(lookup_output_eval_claim)
    }

    #[tracing::instrument(skip_all, name = "LassoNode::verify_claim_reduction")]
    fn verify_claim_reduction(
        &self,
        _: CombinedEvalClaim<E>,
        transcript: &mut dyn TranscriptRead<F, E>,
    ) -> Result<Vec<Vec<EvalClaim<E>>>, Error> {
        let lookup = self.preprocessing.lookups.values().next().unwrap();
        let num_vars = self.num_vars;
        let r = transcript.squeeze_challenges(num_vars);

        let g = self.collation_sum_check_function(lookup, num_vars);
        let claimed_sum = transcript.read_felt_ext()?;

        let _ = info_span!("LassoNode::verify_collation_sum_check")
            .in_scope(|| verify_sum_check(&g, claimed_sum, transcript))?;

        // // Round n+1
        let [gamma, tau] = transcript.squeeze_challenges(2).try_into().unwrap();

        // // memory checking
        Self::verify_memory_checking(&self.preprocessing, num_vars, &gamma, &tau, transcript)?;

        Ok(chain![iter::once(vec![EvalClaim::new(r.to_vec(), claimed_sum)])].collect_vec())
    }
}

impl<F: PrimeField, E: ExtensionField<F>, const C: usize, const M: usize> LassoNode<F, E, C, M> {
    pub fn new(
        // table: Box<dyn DecomposableTable<F, E>>,
        preprocessing: LassoPreprocessing<F, E>,
        num_vars: usize,
        lookups: Vec<LookupId>,
    ) -> Self {
        Self {
            num_vars,
            preprocessing,
            lookups,
        }
    }

    #[tracing::instrument(skip_all, name = "LassoNode::polynomialize")]
    pub fn polynomialize<'a>(
        &self,
        inputs: &BoxMultilinearPoly<F, E>,
    ) -> LassoPolynomials<'a, F, E> {
        let num_reads = inputs.len().next_power_of_two();

        // subtable_lookup_indices : [[usize; num_rows]; num_chunks]
        let subtable_lookup_indices = self.subtable_lookup_indices(inputs);

        let lookup_inputs = izip!(0..inputs.len(), self.lookups.clone())
            .map(|(i, lookup_id)| (i, self.preprocessing.lookup_id_to_index[&lookup_id]))
            .collect_vec();

        let polys: Vec<_> = (0..self.preprocessing.num_memories)
            .into_par_iter()
            .map(|memory_index| {
                let dim_index = self.preprocessing.memory_to_dimension_index[memory_index];
                let subtable_index = self.preprocessing.memory_to_subtable_index[memory_index];
                let access_sequence = &subtable_lookup_indices[dim_index];

                let mut final_cts_i = vec![0usize; M];
                let mut read_cts_i = vec![0usize; num_reads];
                let mut subtable_lookups = vec![F::ZERO; num_reads];

                for (j, lookup) in &lookup_inputs {
                    let memories_used = &self.preprocessing.lookup_to_memory_indices[*lookup];
                    if memories_used.contains(&memory_index) {
                        let memory_address = access_sequence[*j];
                        debug_assert!(memory_address < M);

                        let counter = final_cts_i[memory_address];
                        read_cts_i[*j] = counter;
                        final_cts_i[memory_address] = counter + 1;
                        subtable_lookups[*j] = self
                            .preprocessing
                            .materialized_subtables
                            .as_ref()
                            .expect("subtables not materialized")[subtable_index][memory_address];
                    }
                }

                (
                    DensePolynomial::from_usize::<E>(&read_cts_i),
                    DensePolynomial::from_usize::<E>(&final_cts_i),
                    box_dense_poly(subtable_lookups),
                )
            })
            .collect();

        // Vec<(DensePolynomial<F>, DensePolynomial<F>, DensePolynomial<F>)> -> (Vec<DensePolynomial<F>>, Vec<DensePolynomial<F>>, Vec<DensePolynomial<F>>)
        let (read_cts, final_cts, e_polys) = polys.into_iter().fold(
            (Vec::new(), Vec::new(), Vec::new()),
            |(mut read_acc, mut final_acc, mut e_acc), (read, f, e)| {
                read_acc.push(read);
                final_acc.push(f);
                e_acc.push(e);
                (read_acc, final_acc, e_acc)
            },
        );
        let dims: Vec<_> = subtable_lookup_indices
            .into_par_iter()
            .take(C)
            .map(|mut access_sequence| {
                access_sequence.resize(access_sequence.len().next_power_of_two(), 0);
                DensePolynomial::from_usize(&access_sequence)
            })
            .collect();

        let mut lookup_flag_bitvectors: Vec<Vec<u64>> =
            vec![vec![0u64; num_reads]; self.preprocessing.lookups.len()];

        for (j, lookup_idx) in lookup_inputs.into_iter() {
            lookup_flag_bitvectors[lookup_idx][j] = 1;
        }

        let lookup_flag_polys: Vec<_> = lookup_flag_bitvectors
            .par_iter()
            .map(|flag_bitvector| DensePolynomial::from_u64(flag_bitvector))
            .collect();

        let mut lookup_outputs = self.compute_lookup_outputs(inputs);
        lookup_outputs.resize(num_reads, F::ZERO);
        let lookup_outputs = box_dense_poly(lookup_outputs);

        LassoPolynomials {
            dims,
            read_cts,
            final_cts,
            lookup_flag_polys,
            lookup_flag_bitvectors,
            e_polys,
            lookup_outputs,
        }
    }

    #[allow(clippy::too_many_arguments)]
    #[tracing::instrument(skip_all, name = "LassoNode::prove_collation_sum_check")]
    pub fn prove_collation_sum_check(
        &self,
        lookup_output_poly: &BoxMultilinearPoly<F, E>,
        lookup: &Box<dyn LookupType<F, E>>,
        e_polys: &[BoxMultilinearPoly<F, E>],
        flag_polys: &[BoxMultilinearPoly<F, E>],
        r: &[E],
        num_vars: usize,
        transcript: &mut dyn TranscriptWrite<F, E>,
    ) -> Result<Vec<Vec<EvalClaim<E>>>, Error> {
        let claimed_sum = self.sum_check_claim(r, e_polys, flag_polys);
        if cfg!(feature = "sanity-check") {
            assert_eq!(claimed_sum, lookup_output_poly.evaluate(r));
        }

        transcript.write_felt_ext(&claimed_sum)?;

        let g = self.collation_sum_check_function(lookup, num_vars);

        let polys = e_polys
            .iter()
            .map(|e_poly| SumCheckPoly::Base::<_, _, _, BoxMultilinearPoly<E, E>>(e_poly))
            .collect_vec();

        let (_claim, r_x_prime, e_polys_evals) =
            prove_sum_check(&g, claimed_sum, polys, transcript)?;

        Ok(chain![
            iter::once(vec![EvalClaim::new(r.to_vec(), claimed_sum)]),
            e_polys_evals
                .into_iter()
                .map(|e| vec![EvalClaim::new(r_x_prime.clone(), e)])
        ]
        .collect_vec())
    }

    #[allow(clippy::too_many_arguments)]
    #[tracing::instrument(skip_all, name = "LassoNode::prove_memory_checking")]
    fn prove_memory_checking<'a>(
        preprocessing: &'a LassoPreprocessing<F, E>,
        dims: &'a [BoxMultilinearPoly<'a, F, E>],
        read_ts_polys: &'a [BoxMultilinearPoly<'a, F, E>],
        final_cts_polys: &'a [BoxMultilinearPoly<'a, F, E>],
        e_polys: &'a [BoxMultilinearPoly<'a, F, E>],
        gamma: &E,
        tau: &E,
        transcript: &mut dyn TranscriptWrite<F, E>,
    ) -> Result<(), Error> {
        // key: chunk index, value: chunk
        let mut chunk_map: HashMap<usize, Chunk<F, E>> = HashMap::new();

        let num_memories = preprocessing.num_memories;
        let memories = (0..num_memories).map(|memory_index| {
            let subtable_poly = &preprocessing
                .materialized_subtables
                .as_ref()
                .expect("subtables not materialized")
                [preprocessing.memory_to_subtable_index[memory_index]];
            Memory::<F, E>::new(subtable_poly, &e_polys[memory_index])
        });
        memories.enumerate().for_each(|(memory_index, memory)| {
            let chunk_index = preprocessing.memory_to_dimension_index[memory_index];
            if let std::collections::hash_map::Entry::Vacant(e) = chunk_map.entry(chunk_index) {
                let dim = &dims[chunk_index];
                let read_ts_poly = &read_ts_polys[chunk_index];
                let final_cts_poly = &final_cts_polys[chunk_index];
                e.insert(Chunk::new(
                    chunk_index,
                    dim,
                    read_ts_poly,
                    final_cts_poly,
                    memory,
                ));
            } else {
                chunk_map.entry(chunk_index).and_modify(|chunk| {
                    chunk.add_memory(memory);
                });
            }
        });

        let mut chunks = chunk_map.into_iter().collect_vec();
        chunks.sort_by_key(|(chunk_index, _)| *chunk_index);
        let chunks = chunks.into_iter().map(|(_, chunk)| chunk).collect_vec();

        MemoryCheckingProver::new(chunks, tau, gamma).prove(transcript)
    }

    #[tracing::instrument(skip_all, name = "LassoNode::verify_memory_checking")]
    fn verify_memory_checking(
        preprocessing: &LassoPreprocessing<F, E>,
        num_vars: usize,
        gamma: &E,
        tau: &E,
        transcript: &mut dyn TranscriptRead<F, E>,
    ) -> Result<(), Error> {
        let num_memories = preprocessing.num_memories;
        let mut chunk_map: HashMap<usize, memory_checking::verifier::Chunk<F, E, M>> =
            HashMap::new();
        (0..num_memories).for_each(|memory_index| {
            let chunk_index = preprocessing.memory_to_dimension_index[memory_index];
            let subtable_poly = &preprocessing.subtables_by_idx.as_ref().unwrap()
                [preprocessing.memory_to_subtable_index[memory_index]];
            let memory =
                memory_checking::verifier::Memory::new(memory_index, subtable_poly.clone());

            if let std::collections::hash_map::Entry::Vacant(e) = chunk_map.entry(chunk_index) {
                e.insert(memory_checking::verifier::Chunk::new(
                    chunk_index,
                    log2(M) as usize,
                    memory,
                ));
            } else {
                chunk_map.entry(chunk_index).and_modify(|chunk| {
                    chunk.add_memory(memory);
                });
            }
        });

        let mut chunks = chunk_map.into_iter().collect_vec();
        chunks.sort_by_key(|(chunk_index, _)| *chunk_index);
        let chunks = chunks.into_iter().map(|(_, chunk)| chunk).collect_vec();

        let mem_check = memory_checking::verifier::MemoryCheckingVerifier::new(chunks);
        mem_check.verify(num_vars, gamma, tau, transcript)
    }

    #[tracing::instrument(skip_all, name = "LassoNode::subtable_lookup_indices")]
    fn subtable_lookup_indices(&self, inputs: &BoxMultilinearPoly<F, E>) -> Vec<Vec<usize>> {
        let num_rows: usize = inputs.len();
        let num_chunks = C;

        let indices: Vec<_> = izip_par!((0..num_rows), &self.lookups)
            .map(|(i, lookup_id)| {
                let lookup = &self.preprocessing.lookups[lookup_id];
                let mut index_bits = fe_to_bits_le(inputs[i]);
                index_bits.truncate(lookup.chunk_bits(M).iter().sum());
                if cfg!(feature = "sanity-check") {
                    assert_eq!(
                        usize_from_bits_le(&fe_to_bits_le(inputs[i])),
                        usize_from_bits_le(&index_bits),
                        "index {i} out of range",
                    );
                }
                let mut chunked_index = iter::repeat(0).take(num_chunks).collect_vec();
                let chunked_index_bits = lookup.subtable_indices(index_bits, M.ilog2() as usize);
                chunked_index
                    .iter_mut()
                    .zip(chunked_index_bits)
                    .map(|(chunked_index, index_bits)| {
                        *chunked_index = usize_from_bits_le(&index_bits);
                    })
                    .collect_vec();
                chunked_index
            })
            .collect();

        let lookup_indices = (0..num_chunks)
            .map(|i| indices.iter().map(|indices| indices[i]).collect_vec())
            .collect_vec();
        lookup_indices
    }

    fn compute_lookup_outputs(&self, inputs: &BoxMultilinearPoly<F, E>) -> Vec<F> {
        izip_par!(inputs.as_dense().unwrap(), &self.lookups)
            .map(|(i, lookup_id)| self.preprocessing.lookups[lookup_id].output(i))
            .collect()
    }

    pub fn sum_check_claim(
        &self,
        r: &[E], // claim: CombinedEvalClaim<E>,
        e_polys: &[BoxMultilinearPoly<F, E>],
        flag_polys: &[BoxMultilinearPoly<F, E>],
    ) -> E {
        let num_memories = self.preprocessing.num_memories;
        assert_eq!(e_polys.len(), num_memories);
        let num_vars = e_polys[0].num_vars();
        let bh_size = 1 << num_vars;
        let eq = MultilinearPolynomial::eq_xy(r);
        // \sum_{k \in \{0, 1\}^{\log m}} (\tilde{eq}(r, k) * g(E_1(k), ..., E_{\alpha}(k)))
        let claim = (0..bh_size)
            .into_par_iter()
            .map(|k| {
                eq[k]
                    * izip!(flag_polys, self.preprocessing.lookups.values())
                        .enumerate()
                        .map(|(lookup_idx, (flag_poly, lookup))| {
                            let operands: Vec<_> = self.preprocessing.lookup_to_memory_indices
                                [lookup_idx]
                                .par_iter()
                                .map(|memory_index| e_polys[*memory_index][k])
                                .collect();

                            flag_poly[k] * lookup.combine_lookups(&operands, C, M)
                        })
                        .sum::<F>()
            })
            .sum();

        claim
    }

    // (\tilde{eq}(r, k) * g(E_1(k), ..., E_{\alpha}(k)))
    pub fn collation_sum_check_function(
        &self,
        lookup: &Box<dyn LookupType<F, E>>,
        num_vars: usize,
    ) -> Generic<F, E> {
        let num_memories = self.preprocessing.num_memories;
        // TODO: expression with flag polys
        let exprs = lookup.combine_lookup_expressions(
            (0..num_memories)
                .map(|idx| Expression::poly(idx))
                .collect_vec(),
            C,
            M,
        );

        let eq_r_x = &Expression::poly(0);

        Generic::new(num_vars, &(eq_r_x * exprs))
    }
}

/// All polynomials associated with Jolt instruction lookups.
pub struct LassoPolynomials<'a, F: PrimeField, E: ExtensionField<F>> {
    /// `C` sized vector of `DensePolynomials` whose evaluations correspond to
    /// indices at which the memories will be evaluated. Each `DensePolynomial` has size
    /// `m` (# lookups).
    pub dims: Vec<BoxMultilinearPoly<'a, F, E>>,

    /// `NUM_MEMORIES` sized vector of `DensePolynomials` whose evaluations correspond to
    /// read access counts to the memory. Each `DensePolynomial` has size `m` (# lookups).
    pub read_cts: Vec<BoxMultilinearPoly<'a, F, E>>,

    /// `NUM_MEMORIES` sized vector of `DensePolynomials` whose evaluations correspond to
    /// final access counts to the memory. Each `DensePolynomial` has size M, AKA subtable size.
    pub final_cts: Vec<BoxMultilinearPoly<'a, F, E>>,

    /// `NUM_MEMORIES` sized vector of `DensePolynomials` whose evaluations correspond to
    /// the evaluation of memory accessed at each step of the CPU. Each `DensePolynomial` has
    /// size `m` (# lookups).
    pub e_polys: Vec<BoxMultilinearPoly<'a, F, E>>,

    /// Polynomial encodings for flag polynomials for each instruction.
    /// If using a single instruction this will be empty.
    /// NUM_INSTRUCTIONS sized, each polynomial of length 'm' (# lookups).
    ///
    /// Stored independently for use in sumcheck, combined into single DensePolynomial for commitment.
    pub lookup_flag_polys: Vec<BoxMultilinearPoly<'a, F, E>>,

    /// Instruction flag polynomials as bitvectors, kept in this struct for more efficient
    /// construction of the memory flag polynomials in `read_write_grand_product`.
    pub lookup_flag_bitvectors: Vec<Vec<u64>>,
    /// The lookup output for each instruction of the execution trace.
    pub lookup_outputs: BoxMultilinearPoly<'a, F, E>,
}

#[derive(Debug)]
pub struct LassoPreprocessing<F: Field, E> {
    subtable_to_memory_indices: Vec<Vec<usize>>,
    lookup_to_memory_indices: Vec<Vec<usize>>,
    memory_to_subtable_index: Vec<usize>,
    memory_to_dimension_index: Vec<usize>,
    materialized_subtables: Option<Vec<BoxMultilinearPoly<'static, F, E>>>,
    subtables_by_idx: Option<Vec<Box<dyn LassoSubtable<F, E>>>>,
    lookups: Arc<BTreeMap<LookupId, Box<dyn LookupType<F, E>>>>,
    lookup_id_to_index: HashMap<LookupId, usize>,
    num_memories: usize,
}

impl<F: PrimeField, E: ExtensionField<F>> LassoPreprocessing<F, E> {
    #[tracing::instrument(skip_all, name = "LassoNode::preprocess")]
    pub fn preprocess<const C: usize, const M: usize>(
        lookups: impl IntoIterator<Item = Box<dyn LookupType<F, E>>>,
    ) -> Self {
        let lookups = BTreeMap::from_iter(
            lookups
                .into_iter()
                .map(|lookup| (lookup.lookup_id(), lookup)),
        );

        let lookup_id_to_index: HashMap<_, _> = HashMap::from_iter(
            lookups
                .keys()
                .enumerate()
                .map(|(i, lookup_id)| (lookup_id.clone(), i)),
        );

        let subtables = lookups
            .values()
            .flat_map(|lookup| {
                lookup
                    .subtables(C, M)
                    .into_iter()
                    .map(|(subtable, _)| subtable)
            })
            .unique_by(|subtable| subtable.subtable_id())
            .collect_vec();

        // Build a mapping from subtable type => chunk indices that access that subtable type
        let mut subtable_indices: Vec<SubtableIndices> =
            vec![SubtableIndices::with_capacity(C); subtables.len()];
        let mut subtables_by_idx = vec![None; subtables.len()];
        let mut subtable_id_to_index = HashMap::with_capacity(subtables.len());
        for (_, lookup) in &lookups {
            for (subtable, indices) in lookup.subtables(C, M).into_iter() {
                let subtable_idx = subtable_id_to_index
                    .entry(subtable.subtable_id())
                    .or_insert_with(|| {
                        subtables
                            .iter()
                            .position(|s| s.subtable_id() == subtable.subtable_id())
                            .expect("Subtable not found")
                    });
                subtables_by_idx[*subtable_idx].get_or_insert(subtable);
                subtable_indices[*subtable_idx].union_with(&indices);
            }
        }

        let mut subtable_to_memory_indices = Vec::with_capacity(subtables.len());
        let mut memory_to_subtable_index = vec![];
        let mut memory_to_dimension_index = vec![];

        let mut memory_index = 0;
        for (subtable_index, dimension_indices) in subtable_indices.iter().enumerate() {
            subtable_to_memory_indices
                .push((memory_index..memory_index + dimension_indices.len()).collect_vec());
            memory_to_subtable_index.extend(vec![subtable_index; dimension_indices.len()]);
            memory_to_dimension_index.extend(dimension_indices.iter());
            memory_index += dimension_indices.len();
        }
        let num_memories = memory_index;

        // instruction is a type of lookup
        // assume all instreuctions are the same first
        let mut lookup_to_memory_indices = vec![vec![]; lookups.len()];
        for (lookup_index, lookup_type) in lookups.values().enumerate() {
            for (subtable, dimension_indices) in lookup_type.subtables(C, M) {
                let memory_indices: Vec<_> = subtable_to_memory_indices
                    [subtable_id_to_index[&subtable.subtable_id()]]
                    .iter()
                    .filter(|memory_index| {
                        dimension_indices.contains(memory_to_dimension_index[**memory_index])
                    })
                    .collect();
                lookup_to_memory_indices[lookup_index].extend(memory_indices);
            }
        }

        let materialized_subtables = Some(
            Self::materialize_subtables::<M>(&subtables)
                .into_iter()
                .map(box_dense_poly)
                .collect_vec(),
        );

        Self {
            num_memories,
            materialized_subtables,
            subtable_to_memory_indices,
            memory_to_subtable_index,
            memory_to_dimension_index,
            lookup_to_memory_indices,
            subtables_by_idx: Some(
                subtables_by_idx
                    .into_iter()
                    .map(|s| s.unwrap())
                    .collect_vec(),
            ),
            lookup_id_to_index,
            lookups: Arc::new(lookups),
        }
    }

    fn materialize_subtables<const M: usize>(
        subtables: &[Box<dyn LassoSubtable<F, E>>],
    ) -> Vec<Vec<F>> {
        let mut s = Vec::with_capacity(subtables.len());
        for subtable in subtables.iter() {
            s.push(subtable.materialize(M));
        }
        s
    }

    pub fn to_verifier_preprocessing(&mut self) -> LassoPreprocessing<F, E> {
        LassoPreprocessing {
            subtable_to_memory_indices: self.subtable_to_memory_indices.clone(),
            lookup_to_memory_indices: self.lookup_to_memory_indices.clone(),
            memory_to_subtable_index: self.memory_to_subtable_index.clone(),
            memory_to_dimension_index: self.memory_to_dimension_index.clone(),
            materialized_subtables: None,
            subtables_by_idx: self.subtables_by_idx.take(),
            num_memories: self.num_memories,
            lookup_id_to_index: self.lookup_id_to_index.clone(),
            lookups: self.lookups.clone(),
        }
    }
}
