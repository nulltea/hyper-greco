use super::BfvSkEncryptConstans;

pub struct SkEnc1024_1x27_65537;

impl BfvSkEncryptConstans<1> for SkEnc1024_1x27_65537 {
    const N: usize = 1024;
    const E_BOUND: u64 = 19;
    const S_BOUND: u64 = 1;
    const R1_BOUNDS: [u64; 1] = [1246];
    const R2_BOUNDS: [u64; 1] = [41319090];
    const K1_BOUND: u64 = 32768;
    const QIS: [&str; 1] = ["82638181"];
    const K0IS: [&str; 1] = ["1849798"];
}
