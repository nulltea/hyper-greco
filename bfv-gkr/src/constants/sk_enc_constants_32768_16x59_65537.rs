use super::BfvSkEncryptConstans;

pub struct SkEnc32768_16x59_65537;

impl BfvSkEncryptConstans<16> for SkEnc32768_16x59_65537 {
    const N: usize = 32768;
    const E_BOUND: u64 = 19;
    const S_BOUND: u64 = 1;
    const R1_BOUNDS: [u64; 16] = [17449, 44739, 23337, 23971, 38527, 41231, 47178, 25438, 29323, 21130, 36215, 41226, 47080, 20341, 18613, 40562];
    const R2_BOUNDS: [u64; 16] = [238750987231488128, 238750987231488128, 238750987231488128, 238750987231488128, 238750987231488128, 238750987231488128, 238750987231488128, 238750987231488128, 238750987231488160, 238750987231488160, 238750987231488160, 238750987231488160, 238750987231488160, 238750987231488160, 238750987231488160, 238750987231488192];
    const K1_BOUND: u64 = 32768;
    const QIS: [&'static str; 16] = ["477501974462976263", "477501974462976265", "477501974462976267", "477501974462976269", "477501974462976271", "477501974462976277", "477501974462976283", "477501974462976289", "477501974462976293", "477501974462976299", "477501974462976301", "477501974462976307", "477501974462976311", "477501974462976313", "477501974462976317", "477501974462976361"];
    const K0IS: [&'static str; 16] = ["15519160254606397", "413181248299753132", "101318987089463173", "110550337344198528", "322660099471945682", "362070023329770101", "448729597070750619", "131927434145614106", "188539582116643080", "69158624007851611", "288969678337063080", "361989877431741640", "447294256896967800", "57654044645399264", "32480946673725505", "352321367733214595"];
}