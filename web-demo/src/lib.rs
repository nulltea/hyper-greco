use bfv_gkr::constants;
use bfv_gkr::sk_encryption_circuit::BfvEncrypt;
use goldilocks::Goldilocks;
use goldilocks::GoldilocksExt2;
use rand::thread_rng;
use rayon::prelude::*;
use tracing::info_span;
use wasm_bindgen::prelude::*;

// We need to use the init_thread_pool for it to be publicly visible but it appears unused when
// compiling
#[allow(unused_imports)]
pub use wasm_bindgen_rayon::init_thread_pool;

#[wasm_bindgen]
pub fn init_panic_hook() {
    console_error_panic_hook::set_once();
}

#[wasm_bindgen]
pub struct ProverKey(bfv_gkr::sk_encryption_circuit::ProverKey<Goldilocks, GoldilocksExt2>);

use paste::paste;
use tracing::Level;

#[wasm_bindgen]
pub fn porve_encrypt() {
    let params = BfvParameters::new_with_primes(
        vec![1032193, 1073692673],
        vec![995329, 1073668097],
        40961,
        1 << 11,
    );
    let bounds = bfv_rs::witness_bounds(&params).unwrap();
    let bfv = BfvEncrypt::new(params.clone(), bounds, 1);

    let args = {
        let sk = SecretKey::random_with_params(&params, &mut rng);

        let m: Vec<_> = (0..(params.degree as u64)).collect_vec(); // m here is from lowest degree to largest as input into fhe.rs (REQUIRED)
        let pt = Plaintext::encode(
            &m,
            &params,
            Encoding {
                encoding_type: EncodingType::Poly,
                poly_cache: PolyCache::None,
                level: 0,
            },
        );

        let p = BigInt::from_str_radix("18446744069414584321", 10).unwrap();
        bfv_rs::encrypt_with_witness(params, pt, sk, &mut rng, &p)
            .unwrap()
            .1
    };

    let (pk, vk) = info_span!("setup").in_scope(|| bfv.setup::<Goldilocks, GoldilocksExt2>());
    let proof =
        info_span!("FHE_enc prove").in_scope(|| bfv.prove::<Goldilocks, GoldilocksExt2>(&args, pk));

    let (inputs, _) = info_span!("parse inputs").in_scope(|| bfv.get_inputs(&args));

    info_span!("FHE_enc verify")
        .in_scope(|| bfv.verify::<Goldilocks, GoldilocksExt2>(vk, inputs, args.ct0is, &proof));
}

#[wasm_bindgen]
pub fn parse_args(val: JsValue) -> BfvSkEncryptArgs {
    serde_wasm_bindgen::from_value(val).unwrap()
}

use tracing::level_filters::LevelFilter;
use tracing_forest::tag::NoTag;
use tracing_forest::{ForestLayer, PrettyPrinter, Printer};
use tracing_subscriber::fmt::format::{FmtSpan, Pretty};
use tracing_subscriber::prelude::*;
use tracing_subscriber::EnvFilter;

#[wasm_bindgen(start)]
fn setup() {
    init_panic_hook();

    // For WASM, we must set the directives here at compile time.
    let filter_layer = EnvFilter::default().add_directive(LevelFilter::DEBUG.into());

    let printer = Printer::new().writer(tracing_web::MakeWebConsoleWriter::new());
    let forest_layer = ForestLayer::new(printer, NoTag);
    tracing_subscriber::registry()
        .with(filter_layer)
        .with(forest_layer)
        .init();
    tracing::info!("init!");
}
