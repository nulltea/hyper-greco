use criterion::{
    criterion_group, criterion_main, measurement::Measurement, BenchmarkGroup, BenchmarkId,
    Criterion,
};
use gkr::{
    circuit::node::EvalClaim,
    poly::MultilinearPoly,
    prove_gkr,
    transcript::StdRngTranscript,
    util::{
        arithmetic::{ExtensionField, PrimeField},
        dev::{rand_bytes, rand_vec, seeded_std_rng},
    },
};
use goldilocks::{Goldilocks, GoldilocksExt2};

// fn run_keccak256<F: PrimeField, E: ExtensionField<F>>(
//     field_name: &str,
//     group: &mut BenchmarkGroup<impl Measurement>,
// ) {
//     let setup = |num_reps: usize| {
//         let mut rng = seeded_std_rng();
//         let keccak = Keccak::new(256, num_reps);
//         let input = rand_bytes(num_reps * keccak.rate() - 1, &mut rng);
//         let (_, values) = keccak_circuit(keccak, &input);
//         let output_claims = {
//             let output = &values[74];
//             let point = rand_vec(output.num_vars(), &mut rng);
//             let value = output.evaluate(&point);
//             vec![EvalClaim::new(point, value), EvalClaim::default()]
//         };
//         (input, output_claims)
//     };

//     for num_reps in (5..10).map(|log2| 1 << log2) {
//         let id = BenchmarkId::new(field_name, num_reps);
//         let (input, output_claims) = setup(num_reps);
//         group.bench_with_input(id, &num_reps, |b, _| {
//             b.iter(|| {
//                 let (circuit, values) = keccak_circuit(Keccak::new(256, num_reps), &input);
//                 let mut transcript = StdRngTranscript::default();
//                 prove_gkr::<F, E>(&circuit, &values, &output_claims, &mut transcript).unwrap();
//             });
//         });
//     }
// }

// fn bench_keccak256(c: &mut Criterion) {
//     let mut group = c.benchmark_group("keccak256");
//     group.sample_size(10);

//     run_keccak256::<Goldilocks, GoldilocksExt2>("goldilocks_qe", &mut group);
//     run_keccak256::<bn256::Fr, bn256::Fr>("bn254", &mut group);
// }

// criterion_group!(bench, bench_keccak256);
// criterion_main!(bench);
