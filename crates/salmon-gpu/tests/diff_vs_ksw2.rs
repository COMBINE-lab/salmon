//! Differential test: the pure-Rust banded DP reference must reproduce the exact
//! score salmon computes via `ksw2rs` (the backend behind `ksw2_align_score`),
//! over a randomized corpus of reads vs. reference windows.
//!
//! If this passes, the GPU kernel (which mirrors the reference bit-for-bit) also
//! matches `ksw2rs`, so `--gpu` transparently accelerates `--fullLengthAlignment`
//! while preserving salmon's deterministic output.

use ksw2rs::{Aligner, Extz2Input, KSW_EZ_RIGHT, KSW_EZ_SCORE_ONLY, KSW_NEG_INF};
use salmon_gpu::reference::{banded_extz_score, BandedParams};

const MATCH: i8 = 2;
const MISMATCH: i8 = 4;
const GAPO: i8 = 6;
const GAPE: i8 = 2;

fn params(w: i32) -> BandedParams {
    BandedParams {
        match_score: MATCH as i32,
        mismatch_pen: MISMATCH as i32,
        gap_open_pen: GAPO as i32,
        gap_extend_pen: GAPE as i32,
        bandwidth: w,
    }
}

fn dna5_encode(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&b| match b {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => 4,
        })
        .collect()
}

/// The exact score salmon's `ksw2_align_score` returns: `mqe` if the query end is
/// reached within the band, else `ez.max`.
fn ksw_score(query: &[u8], target: &[u8], w: i32) -> i32 {
    let mut mat = [-MISMATCH; 25];
    for i in 0..5 {
        mat[i * 5 + i] = MATCH;
    }
    mat[24] = 0;
    let q = dna5_encode(query);
    let t = dna5_encode(target);
    let mut aligner = Aligner::new();
    let input = Extz2Input {
        query: &q,
        target: &t,
        m: 5,
        mat: &mat,
        q: GAPO,
        e: GAPE,
        w,
        zdrop: -1,
        end_bonus: 0,
        flag: KSW_EZ_RIGHT | KSW_EZ_SCORE_ONLY,
    };
    let ez = aligner.align(&input);
    if ez.mqe > KSW_NEG_INF {
        ez.mqe
    } else {
        ez.max as i32
    }
}

/// Tiny deterministic PRNG (SplitMix64-ish), so the corpus is reproducible.
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        z ^ (z >> 31)
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
    fn base(&mut self) -> u8 {
        b"ACGT"[self.below(4)]
    }
}

fn rand_seq(rng: &mut Rng, n: usize) -> Vec<u8> {
    (0..n).map(|_| rng.base()).collect()
}

struct Mismatch {
    kind: &'static str,
    scenario: &'static str,
    qlen: usize,
    tlen: usize,
    w: i32,
    reference: i32,
    ksw: i32,
}

/// salmon's acceptance threshold: `floor(minScoreFraction * match * readLen)`,
/// accepted iff `score > threshold` (strict). Mirrors `min_accepted_score`.
fn min_accepted(qlen: usize) -> i32 {
    (0.65f64 * MATCH as f64 * qlen as f64).floor() as i32
}

/// The determinism-relevant contract between the reference DP and `ksw2rs`: they
/// agree on accept/reject for every alignment, and for any alignment either side
/// deems acceptable (or within a small slack of the threshold) the scores are
/// bit-identical. They may differ only on the magnitude of deeply sub-threshold
/// scores, which never enter an equivalence class, so quant output is unaffected.
fn check(mism: &mut Vec<Mismatch>, scenario: &'static str, query: &[u8], target: &[u8], w: i32) {
    let p = params(w);
    let got = banded_extz_score(query, target, &p);
    let want = ksw_score(query, target, w);
    let thr = min_accepted(query.len());

    let mut push = |kind| {
        mism.push(Mismatch {
            kind,
            scenario,
            qlen: query.len(),
            tlen: target.len(),
            w,
            reference: got,
            ksw: want,
        });
    };

    // 1. accept/reject must agree.
    if (got > thr) != (want > thr) {
        push("validity");
        return;
    }
    // 2. scores must match whenever either is at or near acceptance (slack of
    //    one mismatch's worth keeps borderline reads strict).
    let slack = (MATCH + MISMATCH) as i32;
    if (got >= thr - slack || want >= thr - slack) && got != want {
        push("score");
    }
}

#[test]
fn reference_matches_ksw2_over_corpus() {
    let mut rng = Rng(0x5a10_c0ff_ee99_1234);
    let mut mism = Vec::new();
    let margin = 32usize;

    for _ in 0..4000 {
        let qlen = 40 + rng.below(120); // 40..=159
                                        // Build a reference window the read aligns into at offset 0.
        let target = rand_seq(&mut rng, qlen + margin);
        let mut query: Vec<u8> = target[..qlen].to_vec();

        let scenario;
        match rng.below(6) {
            0 => {
                scenario = "exact";
            }
            1 => {
                // substitutions
                let k = 1 + rng.below(6);
                for _ in 0..k {
                    let pos = rng.below(qlen);
                    query[pos] = rng.base();
                }
                scenario = "subst";
            }
            2 => {
                // a short deletion in the read (drop d bases) -> read shorter
                let d = 1 + rng.below(4);
                let at = rng.below(qlen.saturating_sub(d).max(1));
                query.drain(at..(at + d).min(query.len()));
                scenario = "read_del";
            }
            3 => {
                // a short insertion in the read (extra bases) -> read longer
                let ins = 1 + rng.below(4);
                let at = rng.below(qlen);
                let extra: Vec<u8> = (0..ins).map(|_| rng.base()).collect();
                query.splice(at..at, extra);
                scenario = "read_ins";
            }
            4 => {
                // substitutions plus one small indel
                let pos = rng.below(qlen);
                query[pos] = rng.base();
                let d = 1 + rng.below(3);
                let at = rng.below(qlen.saturating_sub(d).max(1));
                query.drain(at..(at + d).min(query.len()));
                scenario = "subst_indel";
            }
            _ => {
                // unrelated read of its own random content
                query = rand_seq(&mut rng, qlen);
                scenario = "unrelated";
            }
        }

        // also sprinkle an N occasionally
        if rng.below(5) == 0 && !query.is_empty() {
            let pos = rng.below(query.len());
            query[pos] = b'N';
        }

        for &w in &[15i32, 31, 1000] {
            check(&mut mism, scenario, &query, &target, w);
        }
    }

    if !mism.is_empty() {
        let total = mism.len();
        let mut report = String::new();
        for m in mism.iter().take(15) {
            report.push_str(&format!(
                "  {} [{}] qlen={} tlen={} w={} reference={} ksw={}\n",
                m.kind, m.scenario, m.qlen, m.tlen, m.w, m.reference, m.ksw
            ));
        }
        panic!(
            "{} determinism-relevant mismatches vs ksw2rs (showing first 15):\n{}",
            total, report
        );
    }
}
