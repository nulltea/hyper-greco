use super::BfvSkEncryptConstans;

pub struct SkEnc2048_1x52_65537;

impl BfvSkEncryptConstans<1> for SkEnc2048_1x52_65537 {
    const N: usize = 2048;
    const E_BOUND: u64 = 19;
    const S_BOUND: u64 = 1;
    const R1_BOUNDS: [u64; 1] = [29882];
    const R2_BOUNDS: [u64; 1] = [1434875766359798];
    const K1_BOUND: u64 = 32768;
    const QIS: [&'static str; 1] = ["2869751532719597"];
    const K0IS: [&'static str; 1] = ["2527239722765942"];
}
