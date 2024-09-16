use super::BfvSkEncryptConstans;

pub struct SkEnc8192_4x55_65537;

impl BfvSkEncryptConstans<4> for SkEnc8192_4x55_65537 {
    const N: usize = 8192;
    const E_BOUND: u64 = 19;
    const S_BOUND: u64 = 1;
    const R1_BOUNDS: [u64; 4] = [28014, 21551, 24676, 29416];
    const R2_BOUNDS: [u64; 4] = [
        13712101976447600,
        13712101976447600,
        13712101976447602,
        13712101976447604,
    ];
    const K1_BOUND: u64 = 32768;
    const QIS: [&'static str; 4] = [
        "27424203952895201",
        "27424203952895203",
        "27424203952895205",
        "27424203952895207",
    ];
    const K0IS: [&'static str; 4] = [
        "20017153978526555",
        "14608220699689817",
        "17223556688605927",
        "21190079862835656",
    ];
}
