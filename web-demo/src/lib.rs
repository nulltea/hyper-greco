use bfv_gkr::constants;
use goldilocks::Goldilocks;
use goldilocks::GoldilocksExt2;
use rand::thread_rng;
use rayon::prelude::*;
use wasm_bindgen::prelude::*;
use bfv_gkr::sk_encryption_circuit::BfvEncrypt;

// We need to use the init_thread_pool for it to be publicly visible but it appears unused when
// compiling
#[allow(unused_imports)]
pub use wasm_bindgen_rayon::init_thread_pool;

#[wasm_bindgen]
pub fn init_panic_hook() {
    console_error_panic_hook::set_once();
}

#[wasm_bindgen]
pub fn parallel_sum(data: &[f64]) -> f64 {
    // Sum the data in parallel
    data.par_iter().sum()
}

#[wasm_bindgen]
pub fn witness_preprocess(n_log2: usize) {
    bfv_rs::gen_witness(n_log2);
}

pub use bfv_gkr::sk_encryption_circuit::BfvSkEncryptArgs;

#[wasm_bindgen]
pub struct ProverKey(bfv_gkr::sk_encryption_circuit::ProverKey<Goldilocks, GoldilocksExt2>);

use paste::paste;

macro_rules! define_bfv_encrypt {
    ($n:expr, $k:expr, $const_type:ty) => {
        paste! {
            #[wasm_bindgen]
            pub struct [<BfvEncrypt $n>](BfvEncrypt<$const_type, $k>);

            #[wasm_bindgen]
            impl [<BfvEncrypt $n>] {
                #[wasm_bindgen]
                pub fn new() -> [<BfvEncrypt $n>] {
                    Self(BfvEncrypt::<$const_type, $k>::new($k))
                }

                #[wasm_bindgen]
                pub fn setup(&self) -> ProverKey {
                    ProverKey(self.0.setup::<Goldilocks, GoldilocksExt2>().0)
                }

                #[wasm_bindgen]
                pub fn prove(
                    &self,
                    args: &BfvSkEncryptArgs,
                    pk: ProverKey,
                ) -> Vec<u8> {
                    self.0.prove::<Goldilocks, GoldilocksExt2>(args, pk.0)
                }
            }
        }
    };
}


define_bfv_encrypt!(1024, 1, constants::SkEnc1024_1x27_65537);
define_bfv_encrypt!(2048, 1, constants::SkEnc2048_1x52_65537);
define_bfv_encrypt!(4096, 2, constants::SkEnc4096_2x55_65537);
define_bfv_encrypt!(8192, 4, constants::SkEnc8192_4x55_65537);
define_bfv_encrypt!(16384, 8, constants::SkEnc16384_8x54_65537);


#[wasm_bindgen]
pub fn parse_args(val: JsValue) -> BfvSkEncryptArgs {
    serde_wasm_bindgen::from_value(val).unwrap()
}