//! Times the alignment kernel in isolation.
//!
//! Alignment is the expensive stage of mapping, but in a real run it is
//! entangled with seeding and chaining, so measuring it alone is the only way to
//! attribute a change to the kernel itself. Both paths salmon can take are
//! measured: the default anchored ksw2 gap/flank cost, and the full-window DP
//! that `--fullLengthAlignment` selects.
//!
//! This replaces a hand-rolled timing loop (`examples/bench_align`). That loop
//! reported one wall-clock average per run, which drifted by up to 28% between
//! runs of the same binary on the same input -- wide enough to hide every
//! few-percent change the mapper has been getting. divan runs each case to a
//! time budget and reports the distribution, so a real change is separable from
//! scheduler noise.
//!
//! Run with `cargo bench -p salmon-map`, or one case with
//! `cargo bench -p salmon-map -- anchored`.

use salmon_map::{align_chain, AlignConfig, Mem, MemChain};

fn main() {
    divan::main();
}

/// Deterministic pseudo-random bases, so every run aligns the same sequences.
fn gen(n: usize, seed: u64) -> Vec<u8> {
    const B: [u8; 4] = *b"ACGT";
    let mut x = seed;
    (0..n)
        .map(|_| {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            B[((x >> 33) & 3) as usize]
        })
        .collect()
}

/// `ncases` (read, reference window, chain) triples: a read taken from the
/// reference with `nmis` substitutions scattered through it, and the two-anchor
/// chain the seeder would have produced for it.
fn batch(read_len: usize, ncases: usize, nmis: usize) -> Vec<(Vec<u8>, Vec<u8>, MemChain)> {
    let mut v = Vec::new();
    for i in 0..ncases {
        let refseq = gen(read_len + 64, 1000 + i as u64);
        let mut read = refseq[0..read_len].to_vec();
        let mut x = 7 + i as u64;
        for _ in 0..nmis {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1);
            let p = (x >> 33) as usize % read_len;
            read[p] = match read[p] {
                b'A' => b'C',
                b'C' => b'G',
                b'G' => b'T',
                _ => b'A',
            };
        }
        let chain = MemChain::new(
            vec![
                Mem::new(0, 0, 25),
                Mem::new((read_len / 2) as i32, (read_len / 2) as i32, 20),
            ],
            45.0,
            true,
        );
        v.push((read, refseq, chain));
    }
    v
}

const READ_LENS: &[usize] = &[51, 100, 150];
const MISMATCHES: &[usize] = &[0, 2];

/// Each sample aligns one batch of 256 cases, so the reported time is per batch;
/// divide by 256 for ns/alignment.
const CASES: usize = 256;

fn run(cases: &[(Vec<u8>, Vec<u8>, MemChain)], cfg: &AlignConfig) -> i64 {
    let mut acc = 0i64;
    for (rd, rf, ch) in cases {
        if let Some(a) = align_chain(rd, rf, ch, cfg) {
            acc += a.score as i64;
        }
    }
    acc
}

#[divan::bench(args = READ_LENS, consts = MISMATCHES)]
fn anchored<const NMIS: usize>(bencher: divan::Bencher, read_len: usize) {
    let cases = batch(read_len, CASES, NMIS);
    let cfg = AlignConfig {
        full_length_alignment: false,
        ..Default::default()
    };
    bencher
        .counter(divan::counter::ItemsCount::new(CASES))
        .bench_local(|| run(divan::black_box(&cases), divan::black_box(&cfg)));
}

#[divan::bench(args = READ_LENS, consts = MISMATCHES)]
fn full_window<const NMIS: usize>(bencher: divan::Bencher, read_len: usize) {
    let cases = batch(read_len, CASES, NMIS);
    let cfg = AlignConfig {
        full_length_alignment: true,
        ..Default::default()
    };
    bencher
        .counter(divan::counter::ItemsCount::new(CASES))
        .bench_local(|| run(divan::black_box(&cases), divan::black_box(&cfg)));
}
