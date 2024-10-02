#[macro_export]
macro_rules! generate_sk_enc_test {
    ($prime_name:expr, $F:ty, $E:ty, $Pcs:ty, $N:expr, $K:expr, $K_BITSIZE:expr) => {
        paste! {
            #[test]
            #[serial_test::serial]
            pub fn [<test_sk_enc_valid_ $prime_name _ $N _ $K x $K_BITSIZE _65537>]() {
                type Params = $crate::constants:: [<SkEnc $N _ $K x $K_BITSIZE _65537 >];
                let env_filter = EnvFilter::builder()
                    .with_default_directive(tracing::Level::INFO.into())
                    .from_env_lossy();

                let subscriber = Registry::default()
                    .with(env_filter)
                    .with(ForestLayer::default());

                let _ = tracing::subscriber::set_global_default(subscriber);

                let rng = seeded_std_rng();

                let file_path = format!(
                    "src/data/{}/sk_enc_{}_{}x{}_65537.json",
                    $prime_name,
                    Params::N,
                    $K,
                    $K_BITSIZE,
                );
                let mut file = File::open(&file_path).expect("Failed to open file");
                let mut data = String::new();
                file.read_to_string(&mut data).expect("Failed to read file");
                let bfv = BfvEncrypt::<Params, $K>::new($K);
                let args =
                    serde_json::from_str::<BfvSkEncryptArgs>(&data).expect("Failed to parse JSON");

                let (pk, vk) =
                    info_span!("setup").in_scope(|| bfv.setup::<$F, $E, $Pcs>(rng.clone()));
                let proof = info_span!("FHE_enc prove")
                    .in_scope(|| bfv.prove::<$F, $E, $Pcs>(&args, pk));

                let (inputs, _) =
                    info_span!("parse inputs").in_scope(|| bfv.get_inputs(&args));

                info_span!("FHE_enc verify")
                    .in_scope(|| bfv.verify::<$F, $E, $Pcs>(vk, inputs, args.ct0is,  &proof));
            }
        }
    };
}
