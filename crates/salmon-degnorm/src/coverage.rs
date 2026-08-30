//! Per-transcript fragment-coverage profiles, accumulated during the mapping
//! pass and written to `aux_info/coverage.gz`.
//!
//! # What is being recorded
//!
//! A coverage profile is "how deep is this transcript covered, as a function of
//! position". Quantification itself never needs it: the EM only cares *how
//! many* fragments a transcript got, not *where* they landed. Degradation does
//! care. RNA degrades from the 5' end inwards, so a degraded molecule is
//! sequenced as a pile of fragments crowded towards its 3' end, and the
//! signature of that is a lopsided coverage curve — the same total count, a
//! different shape.
//!
//! [`CoverageAccumulator`] records that shape while the reads are already being
//! streamed, which is the only cheap moment to do it: the fragment's placement
//! is in hand and the read will never be visited again.
//!
//! # Binning, and why not per base
//!
//! Coverage is accumulated into a fixed number of equal-width bins per
//! transcript (default [`DEFAULT_NUM_BINS`]), not per base. Per-base coverage
//! over a human transcriptome is ~350 M positions, which at 8 bytes is 2.8 GB of
//! live accumulator; 100 bins per transcript is ~160 MB for the same index. The
//! downstream model (see [`crate::nmf`]) fits a smooth envelope over the curve
//! and integrates it, so sub-bin detail would be discarded anyway.
//!
//! Bins are a fixed *fraction* of each transcript, not a fixed number of bases.
//! That is what makes profiles comparable across samples without any further
//! interpolation: bin `b` of transcript `t` means the same stretch of `t` in
//! every sample built against the same index.
//!
//! # The difference-array trick
//!
//! A fragment covers a *range* of bins, so the obvious implementation adds its
//! weight to every bin it touches. On a 2 kb transcript with 100 bins a 200 bp
//! fragment touches ~10 of them, and at tens of millions of fragments × several
//! candidate mappings each that is billions of atomic updates.
//!
//! Instead the accumulator stores the *difference* sequence: adding `w` at the
//! first covered bin and `-w` after the last one, then taking a running sum at
//! the end, reconstructs the same per-bin totals from four updates per mapping
//! regardless of how many bins the fragment spans. The two extra updates (over
//! the naive two) pay for the partially-covered bins at each end, so the result
//! is exact rather than rounded to whole bins.
//!
//! # Fixed point, and why not `AtomicF64`
//!
//! The slots are `AtomicI64` holding a fixed-point value (`FRACTION_BITS`
//! fractional bits) rather than [`salmon_core::AtomicF64`]. `AtomicF64` updates
//! through a compare-and-swap *loop*: under the contention a highly-expressed
//! transcript produces — every worker hitting the same handful of bins — those
//! loops spin. `fetch_add` on an integer is a single uncontended instruction.
//!
//! The cost is quantisation: each update rounds to 2^-20 of a fragment. Weights
//! below [`MIN_MASS`] are dropped before they get there (they are posterior
//! weights of mappings the EM will effectively ignore anyway), and the
//! accumulated bias over even 10^9 fragments stays ~10^-6 relative.

use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;
use std::sync::atomic::{AtomicI64, Ordering};

use anyhow::{bail, Context, Result};

/// Default bins per transcript. Coarse on purpose — see the module docs.
pub const DEFAULT_NUM_BINS: usize = 100;

/// Mappings whose posterior weight is below this contribute nothing. They are
/// also below the fixed-point resolution once split across bins, so this is a
/// clarification of what already happens rather than a new approximation.
pub const MIN_MASS: f64 = 1e-6;

/// Fractional bits of the fixed-point accumulator slots.
const FRACTION_BITS: u32 = 20;
const FIXED_SCALE: f64 = (1u64 << FRACTION_BITS) as f64;

/// File magic + version for the on-disk profile dump.
const MAGIC: &[u8; 8] = b"SLMNCOV1";

/// Lock-free accumulator for per-transcript coverage, shared by every mapping
/// worker.
///
/// Sized once from the index (one row per transcript) and updated with
/// [`add_fragment`](Self::add_fragment) as fragments are placed. Rows are
/// `num_bins + 1` wide: the extra slot absorbs the closing `-w` of a fragment
/// that runs to the very end of the transcript, so the hot path needs no bounds
/// special-case.
#[derive(Debug)]
pub struct CoverageAccumulator {
    num_bins: usize,
    /// `num_bins + 1`; see the struct docs.
    stride: usize,
    ref_lens: Vec<u32>,
    diff: Vec<AtomicI64>,
}

impl CoverageAccumulator {
    /// Allocate an accumulator for transcripts of the given lengths.
    ///
    /// Allocation is `ref_lens.len() * (num_bins + 1) * 8` bytes and is touched
    /// from every worker, so this is only ever built when coverage output was
    /// actually requested.
    pub fn new(ref_lens: Vec<u32>, num_bins: usize) -> Self {
        let num_bins = num_bins.max(1);
        let stride = num_bins + 1;
        let diff = (0..ref_lens.len() * stride)
            .map(|_| AtomicI64::new(0))
            .collect();
        Self {
            num_bins,
            stride,
            ref_lens,
            diff,
        }
    }

    /// Bins per transcript.
    pub fn num_bins(&self) -> usize {
        self.num_bins
    }

    /// Transcripts covered by this accumulator.
    pub fn num_targets(&self) -> usize {
        self.ref_lens.len()
    }

    /// Bytes of accumulator state for `num_targets` transcripts at `num_bins`
    /// bins. Exposed so the CLI can warn before allocating rather than after.
    pub fn footprint_bytes(num_targets: usize, num_bins: usize) -> usize {
        num_targets * (num_bins + 1) * std::mem::size_of::<i64>()
    }

    /// Record one fragment placement: `mass` units of a fragment spanning
    /// `[start, start + span)` on transcript `tid`.
    ///
    /// `mass` is the fragment's posterior weight for this mapping, so an
    /// unambiguous fragment deposits 1.0 and an ambiguous one splits itself
    /// across the transcripts it could have come from. Out-of-range spans are
    /// clipped to the transcript rather than rejected: a proper pair whose
    /// inferred length overruns the annotated end is a real, common event and
    /// the covered prefix is still information.
    pub fn add_fragment(&self, tid: usize, start: i32, span: i32, mass: f64) {
        if mass < MIN_MASS || span <= 0 {
            return;
        }
        let Some(&ref_len) = self.ref_lens.get(tid) else {
            return;
        };
        if ref_len == 0 {
            return;
        }
        let ref_len = f64::from(ref_len);
        let s = f64::from(start.max(0));
        let e = f64::from(start.saturating_add(span)).min(ref_len);
        if e <= s {
            return;
        }

        // Positions expressed in bin units: bin `b` spans `[b, b+1)`.
        let bin_width = ref_len / self.num_bins as f64;
        let x0 = s / bin_width;
        let x1 = e / bin_width;
        let last = self.num_bins - 1;
        let b0 = (x0.floor() as usize).min(last);
        // `ceil() - 1` is the index of the last bin the span touches; a span
        // ending exactly on a bin boundary must not claim the next bin.
        let b1 = ((x1.ceil() as usize).max(1) - 1).min(last);
        let base = tid * self.stride;

        if b0 == b1 {
            // Entirely inside one bin: it is covered over a `x1 - x0` fraction
            // of its width.
            let v = mass * (x1 - x0);
            self.add(base + b0, v);
            self.add(base + b0 + 1, -v);
            return;
        }
        // Partially covered first bin, fully covered middle, partially covered
        // last. Written as differences, the running sum reproduces exactly
        // `v0, 1, 1, …, 1, v1` scaled by `mass`.
        let v0 = (b0 + 1) as f64 - x0;
        let v1 = x1 - b1 as f64;
        self.add(base + b0, mass * v0);
        self.add(base + b0 + 1, mass * (1.0 - v0));
        self.add(base + b1, mass * (v1 - 1.0));
        self.add(base + b1 + 1, -mass * v1);
    }

    fn add(&self, idx: usize, v: f64) {
        let q = (v * FIXED_SCALE).round() as i64;
        if q != 0 {
            self.diff[idx].fetch_add(q, Ordering::Relaxed);
        }
    }

    /// Collapse the difference array into per-bin mean coverage depth.
    ///
    /// A bin's value is the expected number of fragments covering an average
    /// base of that bin: a fragment spanning a whole bin contributes 1.0, one
    /// covering half of it contributes 0.5. Negative values cannot arise from
    /// well-formed input, but fixed-point rounding can leave a slot a few
    /// millionths below zero, so the result is clamped.
    pub fn finish(&self, names: Vec<String>, num_mapped: u64) -> CoverageProfiles {
        assert_eq!(
            names.len(),
            self.ref_lens.len(),
            "coverage profile names must match the accumulator's transcript count"
        );
        let mut values = Vec::with_capacity(self.ref_lens.len() * self.num_bins);
        for t in 0..self.ref_lens.len() {
            let base = t * self.stride;
            let mut running = 0i64;
            for b in 0..self.num_bins {
                running += self.diff[base + b].load(Ordering::Relaxed);
                values.push((running as f64 / FIXED_SCALE).max(0.0) as f32);
            }
        }
        CoverageProfiles {
            num_bins: self.num_bins,
            names,
            ref_lens: self.ref_lens.clone(),
            num_mapped,
            values,
        }
    }
}

/// A finished set of coverage profiles: one row of `num_bins` mean-depth values
/// per transcript, plus the identifying metadata a cross-sample fit needs to
/// check that two samples are talking about the same transcriptome.
#[derive(Debug, Clone, PartialEq)]
pub struct CoverageProfiles {
    pub num_bins: usize,
    pub names: Vec<String>,
    pub ref_lens: Vec<u32>,
    /// Mapped fragments in the run the profiles came from; the sequencing-depth
    /// scale factor for cross-sample comparison.
    pub num_mapped: u64,
    /// Row-major `names.len() × num_bins` mean coverage depth.
    pub values: Vec<f32>,
}

impl CoverageProfiles {
    /// One transcript's profile.
    pub fn row(&self, t: usize) -> &[f32] {
        &self.values[t * self.num_bins..(t + 1) * self.num_bins]
    }

    pub fn num_targets(&self) -> usize {
        self.names.len()
    }

    /// Write the gzipped little-endian dump read back by [`Self::read`].
    ///
    /// The layout is deliberately flat and self-describing — magic, counts,
    /// names, lengths, then one `f32` per (transcript, bin) — so that a
    /// third-party script can consume it without linking salmon. `f32` halves
    /// the file against the `f64` accumulator and still carries ~7 digits,
    /// which is more than a coverage estimate has.
    pub fn write(&self, path: &Path) -> Result<()> {
        let f =
            std::fs::File::create(path).with_context(|| format!("creating {}", path.display()))?;
        let mut w = BufWriter::new(flate2::write::GzEncoder::new(
            f,
            flate2::Compression::default(),
        ));
        w.write_all(MAGIC)?;
        w.write_all(&(self.names.len() as u64).to_le_bytes())?;
        w.write_all(&(self.num_bins as u64).to_le_bytes())?;
        w.write_all(&self.num_mapped.to_le_bytes())?;
        for (name, &len) in self.names.iter().zip(&self.ref_lens) {
            let bytes = name.as_bytes();
            w.write_all(&(bytes.len() as u32).to_le_bytes())?;
            w.write_all(bytes)?;
            w.write_all(&len.to_le_bytes())?;
        }
        // One buffered pass rather than a write per value: at 100 bins over a
        // human transcriptome this is 20 M values.
        let mut buf = Vec::with_capacity(self.values.len() * 4);
        for v in &self.values {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        w.write_all(&buf)?;
        w.flush()?;
        Ok(())
    }

    /// Read a dump written by [`Self::write`].
    pub fn read(path: &Path) -> Result<Self> {
        let f = std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
        let mut r = BufReader::new(flate2::read::GzDecoder::new(f));
        let mut magic = [0u8; 8];
        r.read_exact(&mut magic)
            .with_context(|| format!("reading {}", path.display()))?;
        if &magic != MAGIC {
            bail!(
                "{} is not a salmon coverage dump (bad magic)",
                path.display()
            );
        }
        let num_targets = read_u64(&mut r)? as usize;
        let num_bins = read_u64(&mut r)? as usize;
        let num_mapped = read_u64(&mut r)?;
        if num_bins == 0 {
            bail!("{} declares zero coverage bins", path.display());
        }
        let mut names = Vec::with_capacity(num_targets);
        let mut ref_lens = Vec::with_capacity(num_targets);
        for _ in 0..num_targets {
            let mut n = [0u8; 4];
            r.read_exact(&mut n)?;
            let mut name = vec![0u8; u32::from_le_bytes(n) as usize];
            r.read_exact(&mut name)?;
            names.push(String::from_utf8(name).context("coverage dump: non-UTF-8 target name")?);
            let mut l = [0u8; 4];
            r.read_exact(&mut l)?;
            ref_lens.push(u32::from_le_bytes(l));
        }
        let mut raw = Vec::new();
        r.read_to_end(&mut raw)?;
        let expected = num_targets * num_bins * 4;
        if raw.len() != expected {
            bail!(
                "{}: expected {expected} bytes of coverage values, found {}",
                path.display(),
                raw.len()
            );
        }
        let values = raw
            .chunks_exact(4)
            .map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]]))
            .collect();
        Ok(Self {
            num_bins,
            names,
            ref_lens,
            num_mapped,
            values,
        })
    }
}

fn read_u64<R: BufRead>(r: &mut R) -> Result<u64> {
    let mut b = [0u8; 8];
    r.read_exact(&mut b)?;
    Ok(u64::from_le_bytes(b))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn finish(acc: &CoverageAccumulator, lens: &[u32]) -> Vec<f32> {
        let names = (0..lens.len()).map(|i| format!("t{i}")).collect();
        acc.finish(names, 1).values
    }

    #[test]
    fn a_fragment_spanning_the_whole_transcript_gives_depth_one_everywhere() {
        let acc = CoverageAccumulator::new(vec![1000], 10);
        acc.add_fragment(0, 0, 1000, 1.0);
        for v in finish(&acc, &[1000]) {
            assert!((v - 1.0).abs() < 1e-4, "expected depth 1, got {v}");
        }
    }

    #[test]
    fn a_fragment_covering_one_bin_leaves_the_others_empty() {
        // Bins are 100 b wide; [300, 400) is exactly bin 3.
        let acc = CoverageAccumulator::new(vec![1000], 10);
        acc.add_fragment(0, 300, 100, 1.0);
        let v = finish(&acc, &[1000]);
        for (b, x) in v.iter().enumerate() {
            let want = if b == 3 { 1.0 } else { 0.0 };
            assert!((x - want).abs() < 1e-4, "bin {b}: want {want}, got {x}");
        }
    }

    #[test]
    fn partially_covered_end_bins_get_their_exact_fraction() {
        // Bins are 100 b wide. [250, 620) covers half of bin 2 ([200,300)), all
        // of bins 3–5, and the first fifth of bin 6 ([600,700)).
        let acc = CoverageAccumulator::new(vec![1000], 10);
        acc.add_fragment(0, 250, 370, 1.0);
        let v = finish(&acc, &[1000]);
        let want = [0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 0.2, 0.0, 0.0, 0.0];
        for (b, (&got, &w)) in v.iter().zip(&want).enumerate() {
            assert!((got - w).abs() < 1e-3, "bin {b}: want {w}, got {got}");
        }
    }

    #[test]
    fn mass_scales_the_deposited_depth() {
        let acc = CoverageAccumulator::new(vec![100], 4);
        acc.add_fragment(0, 0, 100, 0.25);
        for v in finish(&acc, &[100]) {
            assert!((v - 0.25).abs() < 1e-4);
        }
    }

    #[test]
    fn spans_running_past_the_transcript_end_are_clipped() {
        let acc = CoverageAccumulator::new(vec![100], 4);
        acc.add_fragment(0, 50, 500, 1.0);
        let v = finish(&acc, &[100]);
        assert!((v[0] - 0.0).abs() < 1e-4);
        assert!((v[1] - 0.0).abs() < 1e-4);
        assert!((v[2] - 1.0).abs() < 1e-4);
        assert!((v[3] - 1.0).abs() < 1e-4);
    }

    #[test]
    fn negligible_and_empty_placements_are_ignored() {
        let acc = CoverageAccumulator::new(vec![100], 4);
        acc.add_fragment(0, 0, 100, MIN_MASS / 2.0);
        acc.add_fragment(0, 0, 0, 1.0);
        acc.add_fragment(9, 0, 50, 1.0); // out-of-range transcript
        assert!(finish(&acc, &[100]).iter().all(|&v| v == 0.0));
    }

    #[test]
    fn rows_are_independent() {
        let acc = CoverageAccumulator::new(vec![100, 200], 4);
        acc.add_fragment(1, 0, 200, 1.0);
        let v = finish(&acc, &[100, 200]);
        assert!(v[..4].iter().all(|&x| x == 0.0));
        assert!(v[4..].iter().all(|x| (x - 1.0).abs() < 1e-4));
    }

    #[test]
    fn profiles_survive_a_write_read_roundtrip() {
        let acc = CoverageAccumulator::new(vec![100, 250], 8);
        acc.add_fragment(0, 10, 60, 0.5);
        acc.add_fragment(1, 100, 150, 1.0);
        let profiles = acc.finish(vec!["txA".into(), "txB".into()], 12345);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("coverage.gz");
        profiles.write(&path).unwrap();
        assert_eq!(CoverageProfiles::read(&path).unwrap(), profiles);
    }

    #[test]
    fn reading_a_non_coverage_file_is_an_error_not_a_panic() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("nope.gz");
        std::fs::write(&path, b"definitely not gzip").unwrap();
        assert!(CoverageProfiles::read(&path).is_err());
    }
}
