//! Timing of the ksw2 alignment validation in salmon's usage scenario:
//! validate a read against the reference window implied by its chain (score-only).
//! Reports the default anchored (ksw2 gap/flank) cost vs the full-window DP
//! (`--fullLengthAlignment`).
use salmon_map::{align_chain, AlignConfig, Mem, MemChain};
use std::time::Instant;

fn gen(n: usize, seed: u64) -> Vec<u8> {
    const B: [u8; 4] = *b"ACGT";
    let mut x = seed;
    (0..n).map(|_| { x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); B[((x>>33)&3) as usize] }).collect()
}

fn batch(read_len: usize, ncases: usize, nmis: usize) -> Vec<(Vec<u8>, Vec<u8>, MemChain)> {
    let mut v = Vec::new();
    for i in 0..ncases {
        let refseq = gen(read_len + 64, 1000 + i as u64);
        let mut read = refseq[0..read_len].to_vec();
        let mut x = 7 + i as u64;
        for _ in 0..nmis {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1);
            let p = (x >> 33) as usize % read_len;
            read[p] = match read[p] { b'A'=>b'C', b'C'=>b'G', b'G'=>b'T', _=>b'A' };
        }
        let chain = MemChain { mems: vec![Mem::new(0, 0, 25), Mem::new((read_len/2) as i32, (read_len/2) as i32, 20)], score: 45.0, is_fw: true };
        v.push((read, refseq, chain));
    }
    v
}

fn time_it(label: &str, cases: &[(Vec<u8>, Vec<u8>, MemChain)], cfg: &AlignConfig, iters: usize) {
    let mut acc = 0i64;
    for (rd, rf, ch) in cases { if let Some(a)=align_chain(rd, rf, ch, cfg) { acc += a.score as i64; } }
    let t = Instant::now();
    for _ in 0..iters {
        for (rd, rf, ch) in cases {
            if let Some(a) = align_chain(rd, rf, ch, cfg) { acc += std::hint::black_box(a.score) as i64; }
        }
    }
    let el = t.elapsed();
    let n = (iters * cases.len()) as f64;
    println!("  {label:32} {:8.1} ns/aln   ({:.2} M aln/s)   [chk {}]", el.as_nanos() as f64 / n, n / el.as_secs_f64() / 1e6, acc & 0xffff);
}

fn main() {
    for &rl in &[51usize, 100, 150] {
        for &nm in &[0usize, 2] {
            let cases = batch(rl, 256, nm);
            let iters = 4000;
            println!("read_len={rl}  mismatches={nm}:");
            let full = AlignConfig { full_length_alignment: true, ..Default::default() };
            let anch = AlignConfig { full_length_alignment: false, ..Default::default() };
            time_it("ksw2 (anchored, DEFAULT path)", &cases, &anch, iters);
            time_it("ksw2 (full-window DP)", &cases, &full, iters);
        }
    }
}
