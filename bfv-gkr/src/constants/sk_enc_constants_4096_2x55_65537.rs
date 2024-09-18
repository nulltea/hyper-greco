use super::BfvSkEncryptConstans;

pub struct SkEnc4096_2x55_65537;

impl BfvSkEncryptConstans<2> for SkEnc4096_2x55_65537 {
    const N: usize = 4096;
    const E_BOUND: u64 = 19;
    const S_BOUND: u64 = 1;
    const R1_BOUNDS: [u64; 2] = [25966, 19503];
    const R2_BOUNDS: [u64; 2] = [13712101976447600, 13712101976447600];
    const K1_BOUND: u64 = 32768;
    const QIS: [&'static str; 2] = ["27424203952895201", "27424203952895203"];
    const K0IS: [&'static str; 2] = ["20017153978526555", "14608220699689817"];
}
