//! Salmon index construction and loading.
//!
//! Salmon's index is a piscem (SSHash + contig table) index built over the
//! compacted de Bruijn graph / reference tiling of the transcriptome. This
//! crate orchestrates the two-stage build that salmon's C++ `BuildSalmonIndex`
//! performs, but using the Rust components:
//!
//! 1. [`cf1_rs::cf_build`] constructs the cDBG and writes the
//!    `.cf_seg` / `.cf_seq` / `.json` tiling.
//! 2. [`piscem_rs::index::build::build_index`] builds the SSHash dictionary and
//!    contig table from that tiling.
//!
//! A small `info.json` records the build parameters. [`SalmonIndex`] wraps the
//! loaded [`piscem_rs`] [`ReferenceIndex`] and exposes the accessors the
//! mapper and quantifier need.

use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use cf1_rs::{cf_build, CfInput};
use piscem_rs::index::build::{build_index, BuildConfig};
use piscem_rs::index::reference_index::ReferenceIndex;
use serde::{Deserialize, Serialize};
use tracing::{info, warn};

/// Basename (within the index directory) of the cDBG tiling files.
const CDBG_PREFIX: &str = "cdbg";
/// Basename (within the index directory) of the piscem index files.
const INDEX_PREFIX: &str = "index";
/// Name of the salmon index metadata file.
const INFO_FILE: &str = "info.json";
/// Concatenated reference sequences (forward strand), in transcript-id order.
const REFSEQ_FILE: &str = "refseq.bin";
/// Per-reference offsets into [`REFSEQ_FILE`] (`num_refs + 1` entries).
const REFSEQ_OFFSETS_FILE: &str = "refseq_offsets.json";

/// Default minimizer length for a given k (used when `m == 0`).
fn default_minimizer_len(k: usize) -> usize {
    (k / 2) + 1
}

/// Resolve a thread count, mapping `0` to all available cores.
fn effective_threads(threads: usize) -> usize {
    if threads == 0 {
        std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1)
    } else {
        threads
    }
}

/// Options for building a salmon index.
#[derive(Debug, Clone)]
pub struct IndexBuildOptions {
    /// transcriptome FASTA file(s)
    pub transcripts: Vec<PathBuf>,
    /// output index directory (created if absent)
    pub output_dir: PathBuf,
    /// k-mer length (must be odd, ≤ 63); salmon's default is 31
    pub k: usize,
    /// minimizer length; `0` selects [`default_minimizer_len`]
    pub m: usize,
    /// worker threads; `0` means all available cores
    pub threads: usize,
    /// build canonical-k-mer index (the salmon/pufferfish convention)
    pub canonical: bool,
    /// build the equivalence-class table (needed for pseudoalignment-only mode)
    pub build_ec_table: bool,
    /// keep the intermediate cDBG `.cf_*` files instead of deleting them
    pub keep_intermediate: bool,
    /// keep the post-cleaning ("fixed") reference FASTA — the copy with non-ACGT
    /// bases replaced — instead of deleting it (salmon's `--keepFixedFasta`).
    /// Implied by `keep_intermediate`.
    pub keep_fixed_fasta: bool,
    /// directory for build intermediates (the cleaned FASTA + cDBG `.cf_*`
    /// files) instead of the output dir (salmon's `--tmpdir`); the final index
    /// still lands in `output_dir`. Created if absent.
    pub tmpdir: Option<PathBuf>,
    /// retain exact-duplicate transcript sequences instead of collapsing them
    /// (salmon's `--keepDuplicates`; default `false` = remove duplicates).
    pub keep_duplicates: bool,
    /// optional file of decoy sequence names (one per line; salmon's `--decoys`).
    /// Records whose name is in this set are indexed but flagged as decoys; they
    /// must appear contiguously at the end of the FASTA input.
    pub decoys: Option<PathBuf>,
    /// GENCODE-style headers: additionally truncate each transcript name at the
    /// first `|` (salmon's `--gencode`), so e.g. `ENST...|ENSG...|...` becomes
    /// `ENST...`. Applied consistently to the index names and the seq/name hashes.
    pub gencode: bool,
    /// Clip poly-A tails from reference ends before indexing, matching salmon/
    /// pufferfish `FixFasta` (default `true`; salmon's `--no-clip`/`-n` disables).
    /// A reference whose (cleaned) sequence ends in a run of >= 10 `A`s has all
    /// its trailing `A`s trimmed; an all-`A` reference is dropped from the index.
    /// Reference hashes are computed pre-clip, so this does not change provenance.
    pub clip_polya: bool,
}

impl IndexBuildOptions {
    /// Construct options with salmon's defaults for the given transcripts and
    /// output directory.
    pub fn new(transcripts: Vec<PathBuf>, output_dir: PathBuf) -> Self {
        Self {
            transcripts,
            output_dir,
            k: 31,
            m: 0,
            threads: 0,
            canonical: true,
            build_ec_table: true,
            keep_intermediate: false,
            keep_fixed_fasta: false,
            tmpdir: None,
            keep_duplicates: false,
            decoys: None,
            gencode: false,
            clip_polya: true,
        }
    }
}

/// Metadata written to `info.json` alongside the piscem index.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IndexInfo {
    pub k: usize,
    pub m: usize,
    pub canonical: bool,
    pub has_ec_table: bool,
    pub num_refs: usize,
    /// Index of the first decoy reference: references `[first_decoy_index, num_refs)`
    /// are decoys (genome/contamination sequences indexed for selective-alignment
    /// but never quantified). `None` (or absent) means no decoys. salmon requires
    /// all decoy records to appear contiguously at the end of the FASTA.
    #[serde(default)]
    pub first_decoy_index: Option<usize>,
    /// Number of decoy references (the contiguous block beginning at
    /// `first_decoy_index`). Recorded separately because the index builder
    /// (piscem) appends sub-`k` "short" references *after* the decoy block, so
    /// decoys are **not** the suffix of the reference numbering: transcripts
    /// occupy `[0, first_decoy_index)`, decoys `[first_decoy_index,
    /// first_decoy_index + num_decoys)`, and short (never-seeded, 0-count)
    /// transcripts the remaining tail `[first_decoy_index + num_decoys,
    /// num_refs)`. `0` when there are no decoys.
    #[serde(default)]
    pub num_decoys: usize,
    /// salmon version that produced the index
    pub salmon_version: String,
    /// SHA-256/512 of the reference sequences and names, computed to byte-match
    /// salmon's pufferfish `FixFasta` hashes (so downstream provenance tools
    /// recognize a Rust-built index as equivalent for the same FASTA). Empty for
    /// indexes built before these were recorded (`#[serde(default)]`).
    #[serde(default)]
    pub seq_hash: String,
    #[serde(default)]
    pub name_hash: String,
    #[serde(default)]
    pub seq_hash512: String,
    #[serde(default)]
    pub name_hash512: String,
    #[serde(default)]
    pub decoy_seq_hash: String,
    #[serde(default)]
    pub decoy_name_hash: String,
}

/// Compute the reference sequence/name SHA-256 and SHA-512 hashes exactly as
/// salmon's pufferfish `FixFasta` does, so the values byte-match a salmon index
/// built from the same FASTA(s):
/// - sequences are hashed in FASTA order, with non-printable bytes removed but
///   **before** uppercasing / non-ACGT replacement / poly-A clipping (i.e. the
///   original-case printable sequence), for every record;
/// - names are the header up to the first space or tab (salmon's `sepStr=" \t"`),
///   hashed for every record whose (printable) sequence is non-empty;
/// - decoy hashers receive nothing here (decoys are not yet supported), so the
///   decoy digests are the SHA of the empty string, matching salmon's no-decoy case.
fn compute_ref_hashes(
    transcripts: &[PathBuf],
    decoy_names: &ahash::AHashSet<Vec<u8>>,
    gencode: bool,
) -> Result<RefHashes> {
    use sha2::{Digest, Sha256, Sha512};
    let mut seq256 = Sha256::new();
    let mut name256 = Sha256::new();
    let mut seq512 = Sha512::new();
    let mut name512 = Sha512::new();
    // Decoy sequences/names are routed to their own hashers (salmon's
    // `update_seq_hash(isDecoy,..)` / `update_name_hash`); with no decoys these
    // stay empty → SHA of "" (matches salmon's no-decoy case). Only the SHA-256
    // decoy digests are recorded by salmon, so we only need those.
    let mut decoy_seq256 = Sha256::new();
    let mut decoy_name256 = Sha256::new();

    for path in transcripts {
        let mut reader = needletail::parse_fastx_file(path)
            .with_context(|| format!("opening {} for hashing", path.display()))?;
        while let Some(rec) = reader.next() {
            let rec = rec.context("reading FASTA record for hashing")?;
            // Keep only C-`isprint` bytes (0x20..=0x7e); preserve original case.
            let seq: Vec<u8> = rec
                .seq()
                .iter()
                .copied()
                .filter(|&b| (0x20..=0x7e).contains(&b))
                .collect();
            // Header up to the first ' '/'\t' (and first '|' under --gencode).
            let id = rec.id();
            let name = processed_name(id, gencode);
            let is_decoy = decoy_names.contains(name);
            if is_decoy {
                decoy_seq256.update(&seq);
            } else {
                seq256.update(&seq);
                seq512.update(&seq);
            }
            if !seq.is_empty() {
                if is_decoy {
                    decoy_name256.update(name);
                } else {
                    name256.update(name);
                    name512.update(name);
                }
            }
        }
    }
    let hex = |bytes: &[u8]| -> String {
        let mut s = String::with_capacity(bytes.len() * 2);
        for b in bytes {
            s.push_str(&format!("{b:02x}"));
        }
        s
    };
    Ok(RefHashes {
        seq_hash: hex(&seq256.finalize()),
        name_hash: hex(&name256.finalize()),
        seq_hash512: hex(&seq512.finalize()),
        name_hash512: hex(&name512.finalize()),
        decoy_seq_hash: hex(&decoy_seq256.finalize()),
        decoy_name_hash: hex(&decoy_name256.finalize()),
    })
}

struct RefHashes {
    seq_hash: String,
    name_hash: String,
    seq_hash512: String,
    name_hash512: String,
    decoy_seq_hash: String,
    decoy_name_hash: String,
}

/// Stream the input FASTA(s) replacing any base that is not A/C/G/T (case-
/// insensitive) with a pseudo-random ACGT base, writing one cleaned FASTA to
/// `out_path` (record order and headers preserved). salmon's `FixFasta` does the
/// same: cf1-rs/packed-seq cannot index `N`/ambiguous bases. The reference seq/
/// name HASHES are computed on the *original* input (matching salmon, which
/// hashes before replacement). Returns the number of bases replaced.
/// Outcome of [`preprocess_fasta`]: how many bases were randomized and which
/// transcripts were dropped as exact-sequence duplicates.
struct PreprocessResult {
    /// non-ACGT bases replaced with pseudo-random ACGT
    replaced: u64,
    /// `(retained_name, duplicate_name)` pairs in discovery order; the duplicate
    /// was *not* written to the cleaned FASTA (so it is absent from the index).
    duplicate_clusters: Vec<(Vec<u8>, Vec<u8>)>,
    /// number of references whose trailing poly-A run was clipped.
    polya_clipped: u64,
    /// names of references dropped because clipping left them empty (all-`A`).
    polya_dropped: Vec<Vec<u8>>,
    /// names of *decoy* references dropped because their cleaned/clipped length is
    /// `<= k`: such a reference carries no k-mers, can never be seeded, and so can
    /// never function as a decoy. Keeping it would only let the index builder
    /// relocate it into the trailing sub-`k` block, breaking the contiguity of the
    /// decoy region that decoy classification relies on. (Sub-`k` *transcripts*
    /// are kept — they are still reported in `quant.sf` with 0 reads.)
    short_decoys_dropped: Vec<Vec<u8>>,
}

/// salmon/pufferfish `FixFasta` clips a reference iff its (cleaned) sequence ends
/// in a run of at least this many `A`s; the entire trailing `A` run is then
/// removed (an all-`A` reference is dropped).
const POLYA_CLIP_LEN: usize = 10;

/// Length of the longest run of consecutive ACGT (case-insensitive) bases — the
/// longest stretch Cuttlefish can extract k-mers from (it splits the de Bruijn
/// graph on any other base, e.g. `N`). A run of length `L` carries `L - k + 1`
/// k-mers, so a reference is tileable iff its longest ACGT run is `>= k`.
fn longest_acgt_run(seq: &[u8]) -> usize {
    let mut best = 0usize;
    let mut cur = 0usize;
    for &b in seq {
        if matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't') {
            cur += 1;
            best = best.max(cur);
        } else {
            cur = 0;
        }
    }
    best
}

fn preprocess_fasta(
    inputs: &[PathBuf],
    out_path: &Path,
    keep_duplicates: bool,
    decoy_names: &ahash::AHashSet<Vec<u8>>,
    gencode: bool,
    clip_polya: bool,
    k: usize,
) -> Result<PreprocessResult> {
    use std::io::Write as _;
    const B: [u8; 4] = *b"ACGT";
    let mut x: u64 = 0x9E37_79B9_7F4A_7C15;
    let mut out = std::io::BufWriter::new(
        std::fs::File::create(out_path)
            .with_context(|| format!("creating {}", out_path.display()))?,
    );
    let mut replaced = 0u64;
    let mut polya_clipped = 0u64;
    let mut polya_dropped: Vec<Vec<u8>> = Vec::new();
    let mut short_decoys_dropped: Vec<Vec<u8>> = Vec::new();
    let mut duplicate_clusters: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
    let mut saw_decoy = false;
    // Exact original (pure-ACGT) sequence -> name of the first transcript that
    // carried it. salmon collapses transcripts whose *processed* sequences match;
    // because it randomizes non-ACGT bases with a position-advancing RNG, two
    // identical sequences that contain any non-ACGT base never compare equal and
    // so are never collapsed — we mirror that by only deduplicating pure-ACGT
    // sequences (keyed on the exact original bytes, matching salmon's pre-process
    // XXH64 over the unmodified sequence). Only populated when removing dups.
    let mut seen: ahash::AHashMap<Vec<u8>, Vec<u8>> = ahash::AHashMap::new();
    for path in inputs {
        let mut reader = needletail::parse_fastx_file(path)
            .with_context(|| format!("opening {}", path.display()))?;
        while let Some(rec) = reader.next() {
            let rec = rec.context("reading FASTA record")?;
            // Processed name: header up to the first space/tab (salmon's sepStr),
            // and the first `|` too under `--gencode`.
            let id = rec.id();
            let name = processed_name(id, gencode);
            let orig = rec.seq();
            let is_decoy = decoy_names.contains(name);
            // salmon requires every decoy record to appear, contiguously, at the
            // end of the input — reject a real transcript after any decoy.
            if is_decoy {
                saw_decoy = true;
            } else if saw_decoy {
                anyhow::bail!(
                    "non-decoy reference {:?} appears after a decoy; decoy records must be \
                     contiguous at the end of the FASTA input",
                    String::from_utf8_lossy(name)
                );
            }

            if !keep_duplicates {
                let pure_acgt = orig
                    .iter()
                    .all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'));
                if pure_acgt {
                    if let Some(retained_name) = seen.get(orig.as_ref()) {
                        // Exact duplicate of an earlier reference: record the
                        // cluster and drop it from the index entirely.
                        duplicate_clusters.push((retained_name.clone(), name.to_vec()));
                        continue;
                    }
                    seen.insert(orig.to_vec(), name.to_vec());
                }
            }

            let mut seq = orig.into_owned();
            // Replace non-ACGT bases with pseudo-random ACGT (salmon FixFasta
            // behavior). This is currently applied to decoys too: ideally decoys
            // would keep their N-runs so Cuttlefish splits the graph on them
            // (yielding a much smaller/faster cDBG for genome decoys with large
            // assembly gaps), but cf1-rs (<= 0.4.1) panics on `N` in its
            // minimizer-counting phase (packed-seq rejects non-ACGTN) despite
            // `poly_n_stretch` being set — so we cannot yet feed it N. The
            // longest-ACGT-run decoy-drop below is already written to that future
            // (it reduces to a length test while decoys are N-free).
            for b in seq.iter_mut() {
                if !matches!(*b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't') {
                    x = x
                        .wrapping_mul(6364136223846793005)
                        .wrapping_add(1442695040888963407);
                    *b = B[((x >> 33) & 3) as usize];
                    replaced += 1;
                }
            }
            // Kallisto-esque poly-A clipping (salmon/pufferfish FixFasta): if the
            // (cleaned) sequence ends in a run of >= POLYA_CLIP_LEN `A`s, trim all
            // trailing `A`s. This runs after non-ACGT replacement (matching salmon's
            // order); the reference hashes are computed pre-clip, so provenance is
            // unaffected.
            if clip_polya
                && seq.len() > POLYA_CLIP_LEN
                && seq[seq.len() - POLYA_CLIP_LEN..]
                    .iter()
                    .all(|&b| matches!(b, b'A' | b'a'))
            {
                match seq.iter().rposition(|&b| !matches!(b, b'A' | b'a')) {
                    Some(last) => {
                        seq.truncate(last + 1);
                        polya_clipped += 1;
                    }
                    None => {
                        // All `A`s: drop the reference entirely (not written, not
                        // counted as retained), matching salmon.
                        polya_clipped += 1;
                        polya_dropped.push(name.to_vec());
                        continue;
                    }
                }
            }

            // Drop a decoy that Cuttlefish cannot tile. A reference carries a k-mer
            // (and so a unitig) iff it has a run of >= k consecutive ACGT bases — a
            // run of exactly k has one k-mer (validated empirically). A decoy with
            // no such run (too short, or fragmented into sub-k pieces by N-runs)
            // can never be seeded or catch a read, and keeping it would let the
            // builder relocate it into the trailing sub-`k`/untiled block, breaking
            // the contiguity that the O(1) decoy check relies on. Sub-`k`/N-broken
            // *transcripts* are kept (reported at 0 reads); after the replacement
            // above transcripts have no N-runs, so only decoys are affected here.
            if is_decoy && longest_acgt_run(&seq) < k {
                short_decoys_dropped.push(name.to_vec());
                continue;
            }

            out.write_all(b">")?;
            out.write_all(name)?;
            out.write_all(b"\n")?;
            out.write_all(&seq)?;
            out.write_all(b"\n")?;
        }
    }
    out.flush()?;
    Ok(PreprocessResult {
        replaced,
        duplicate_clusters,
        polya_clipped,
        polya_dropped,
        short_decoys_dropped,
    })
}

/// Read a decoy-names file (one name per line; blank lines ignored) into a set,
/// matching how `name` is extracted in [`preprocess_fasta`] (no `>` prefix).
fn read_decoy_names(path: &Path) -> Result<ahash::AHashSet<Vec<u8>>> {
    let bytes =
        std::fs::read(path).with_context(|| format!("reading decoys file {}", path.display()))?;
    let mut set = ahash::AHashSet::new();
    for line in bytes.split(|&b| b == b'\n') {
        // trim trailing CR and surrounding whitespace
        let line: &[u8] = {
            let mut s = line;
            while let [rest @ .., last] = s {
                if matches!(last, b'\r' | b' ' | b'\t') {
                    s = rest;
                } else {
                    break;
                }
            }
            while let [first, rest @ ..] = s {
                if matches!(first, b' ' | b'\t') {
                    s = rest;
                } else {
                    break;
                }
            }
            s
        };
        if !line.is_empty() {
            set.insert(line.to_vec());
        }
    }
    Ok(set)
}

/// Write salmon's `duplicate_clusters.tsv` (`RetainedRef<TAB>DuplicateRef`, one
/// line per collapsed duplicate). Always written (salmon emits it even with a
/// header and no rows) so downstream tools find the expected file.
fn write_duplicate_clusters(dir: &Path, clusters: &[(Vec<u8>, Vec<u8>)]) -> Result<()> {
    use std::io::Write as _;
    let path = dir.join("duplicate_clusters.tsv");
    let mut f = std::io::BufWriter::new(
        std::fs::File::create(&path).with_context(|| format!("creating {}", path.display()))?,
    );
    f.write_all(b"RetainedRef\tDuplicateRef\n")?;
    for (retained, dup) in clusters {
        f.write_all(retained)?;
        f.write_all(b"\t")?;
        f.write_all(dup)?;
        f.write_all(b"\n")?;
    }
    f.flush()?;
    Ok(())
}

/// Build a salmon index from a transcriptome FASTA into `opts.output_dir`.
pub fn build(opts: &IndexBuildOptions) -> Result<IndexInfo> {
    anyhow::ensure!(!opts.transcripts.is_empty(), "no transcript files provided");
    anyhow::ensure!(
        opts.k >= 1 && opts.k <= 63 && opts.k % 2 == 1,
        "k must be odd in [1, 63]"
    );

    std::fs::create_dir_all(&opts.output_dir)
        .with_context(|| format!("creating index dir {}", opts.output_dir.display()))?;

    let threads = effective_threads(opts.threads);
    let m = if opts.m == 0 {
        default_minimizer_len(opts.k)
    } else {
        opts.m
    };
    // Build intermediates (cleaned FASTA + cDBG files) go to `tmpdir` when set,
    // else the output dir; the final index always lands in the output dir.
    let intermediate_dir = match &opts.tmpdir {
        Some(d) => {
            std::fs::create_dir_all(d)
                .with_context(|| format!("creating index tmpdir {}", d.display()))?;
            d.clone()
        }
        None => opts.output_dir.clone(),
    };
    let cdbg_prefix = intermediate_dir.join(CDBG_PREFIX);
    let index_prefix = opts.output_dir.join(INDEX_PREFIX);

    // Decoy names (genome/contamination sequences indexed but never quantified).
    let decoy_names = match &opts.decoys {
        Some(p) => read_decoy_names(p)?,
        None => ahash::AHashSet::new(),
    };

    // Preprocess the references into the cleaned FASTA that feeds both the cDBG
    // build and the refseq store (the hashes below use the original input):
    // transcripts get non-ACGT bases replaced with pseudo-random ACGT (salmon
    // FixFasta behavior); decoys keep their non-ACGT bases (e.g. genome N-runs) so
    // Cuttlefish splits the graph on them instead of tiling random replacements.
    let cleaned = intermediate_dir.join("cleaned_refs.fa");
    let pre = preprocess_fasta(
        &opts.transcripts,
        &cleaned,
        opts.keep_duplicates,
        &decoy_names,
        opts.gencode,
        opts.clip_polya,
        opts.k,
    )
    .context("preprocessing reference FASTA (non-ACGT replacement, dedup)")?;
    if pre.replaced > 0 {
        info!(
            "replaced {} non-ACGT base(s) with random bases",
            pre.replaced
        );
    }
    if pre.polya_clipped > 0 {
        info!(
            "clipped poly-A tails from {} transcripts",
            pre.polya_clipped
        );
    }
    for dropped in &pre.polya_dropped {
        warn!(
            "transcript {:?} was all A's after cleaning; dropped from the index",
            String::from_utf8_lossy(dropped)
        );
    }
    if !pre.short_decoys_dropped.is_empty() {
        warn!(
            "dropped {} decoy reference(s) with no run of >= k={} consecutive ACGT bases \
             (too short, or fragmented by N-runs): they carry no k-mers and can never be \
             seeded or catch a read: {}",
            pre.short_decoys_dropped.len(),
            opts.k,
            pre.short_decoys_dropped
                .iter()
                .map(|n| String::from_utf8_lossy(n).into_owned())
                .collect::<Vec<_>>()
                .join(", ")
        );
    }
    if !opts.keep_duplicates {
        if !pre.duplicate_clusters.is_empty() {
            info!(
                "removed {} transcripts that were exact sequence duplicates \
                 (use --keepDuplicates to retain them)",
                pre.duplicate_clusters.len()
            );
        }
        write_duplicate_clusters(&opts.output_dir, &pre.duplicate_clusters)
            .context("writing duplicate_clusters.tsv")?;
    }

    // Stage 1: compacted de Bruijn graph / reference tiling.
    info!("building compacted dBG (k={}, threads={})", opts.k, threads);
    let cf = cf_build()
        .input(CfInput::Files(vec![cleaned.clone()]))
        .output_prefix(cdbg_prefix.clone())
        .k(opts.k)
        .threads(threads)
        // Emit references in input order so the decoy block (always last in the
        // input) stays contiguous in the built index — letting `is_decoy` be a
        // single O(1) range test and making the build deterministic.
        .synchronize_output(true)
        .call()
        .context("cf1-rs cDBG construction failed")?;
    info!(
        "cDBG: {} unitigs, {} distinct k-mers",
        cf.unitig_count, cf.vertex_count
    );
    // Authoritative cross-check of our pre-detection against Cuttlefish's own
    // tiling: having dropped every decoy with no >= k ACGT run at preprocess,
    // Cuttlefish should report no untiled sequence (transcripts have no N-runs
    // after replacement; sub-k transcripts are "short", not "untiled"). A
    // non-empty list means our longest-ACGT-run criterion disagreed with
    // Cuttlefish for some reference — warn here; the decoy-contiguity guard below
    // is the hard backstop against any resulting mis-classification.
    if !cf.untiled_seqs.is_empty() {
        warn!(
            "Cuttlefish reported {} untiled reference(s) not dropped at preprocess: {:?}. \
             This is unexpected (decoys with no >= k ACGT run should already be dropped); \
             please report.",
            cf.untiled_seqs.len(),
            cf.untiled_seqs
                .iter()
                .map(|(n, _)| n.as_str())
                .collect::<Vec<_>>()
        );
    }

    // Stage 2: piscem SSHash + contig-table index over the tiling.
    info!("building piscem index (m={})", m);
    let config = BuildConfig {
        input_prefix: cdbg_prefix.clone(),
        output_prefix: index_prefix.clone(),
        k: opts.k,
        m,
        build_ec_table: opts.build_ec_table,
        num_threads: opts.threads,
        canonical: opts.canonical,
        seed: 1,
        single_mphf: false,
        emit_tiny: None,
    };
    build_index(&config).context("piscem index construction failed")?;

    // Load once to capture reference count and confirm the index is usable.
    let idx = ReferenceIndex::load(&index_prefix, opts.build_ec_table, false)
        .context("loading freshly built index")?;
    // Reference seq/name hashes (salmon-FixFasta-compatible), over the input FASTA.
    let h = compute_ref_hashes(&opts.transcripts, &decoy_names, opts.gencode)
        .context("hashing reference sequences/names")?;
    // Locate the decoy block in the *built* index's reference numbering. We
    // cannot reuse the preprocess `retained` count (`pre.first_decoy_index`):
    // piscem moves sub-`k` "short" references — which carry no k-mers and so
    // never tile into the de Bruijn graph — to the *end* of the reference list,
    // after the decoy block. The preprocess count includes those shorts among
    // the transcripts, so it overshoots the decoys' true start by the number of
    // shorts. Scan the actual names instead; this keeps `is_decoy`, the
    // bias-model bounds, and decoy-read accounting correct, and prevents a 250 Mb
    // genome decoy from being swept as a "transcript" in the expected-bias build.
    //
    // Decoys form one contiguous block `[first_decoy_index, first_decoy_index +
    // num_decoys)` (the genome sits after all transcripts, before the trailing
    // sub-`k` shorts). Dropping sub-`k` decoys at preprocess is what keeps that
    // block contiguous: a relocated short decoy would otherwise scatter into the
    // trailing block. We rely on contiguity so the per-fragment decoy check is a
    // single O(1) range test (matching C++ salmon's `tid >= firstDecoyIndex`),
    // rather than a per-mapping set lookup. As a build-time safety net we verify
    // contiguity here and fail loudly if it is ever violated, rather than silently
    // mis-classifying a transcript as a decoy (or vice-versa) at quant time.
    let (first_decoy_index, num_decoys) = if decoy_names.is_empty() {
        (None, 0usize)
    } else {
        let n = idx.num_refs();
        let decoy_tids: Vec<usize> = (0..n)
            .filter(|&t| decoy_names.contains(idx.ref_name(t).as_bytes()))
            .collect();
        match decoy_tids.first().copied() {
            None => (None, 0usize),
            Some(first) => {
                let count = decoy_tids.len();
                // contiguous iff the decoy tids are exactly first..first+count
                let contiguous = decoy_tids
                    .iter()
                    .enumerate()
                    .all(|(i, &t)| t == first + i);
                if !contiguous {
                    anyhow::bail!(
                        "decoy references are not contiguous in the built index \
                         (first decoy at {first}, {count} decoys, but they do not occupy \
                         {first}..{}). This is unexpected after sub-k decoys are dropped at \
                         preprocess — please report, including your decoy set. Built-index \
                         decoy indices: {decoy_tids:?}",
                        first + count
                    );
                }
                (Some(first), count)
            }
        }
    };
    if let Some(fdi) = first_decoy_index {
        let n_short = idx.num_refs() - (fdi + num_decoys);
        info!(
            "{num_decoys} decoy reference(s) (first decoy index {fdi}); \
             {n_short} short (sub-k, 0-count) transcript(s) recorded after the decoys"
        );
    }
    let info = IndexInfo {
        k: opts.k,
        m,
        canonical: opts.canonical,
        has_ec_table: idx.has_ec_table(),
        num_refs: idx.num_refs(),
        first_decoy_index,
        num_decoys,
        salmon_version: env!("CARGO_PKG_VERSION").to_string(),
        seq_hash: h.seq_hash,
        name_hash: h.name_hash,
        seq_hash512: h.seq_hash512,
        name_hash512: h.name_hash512,
        decoy_seq_hash: h.decoy_seq_hash,
        decoy_name_hash: h.decoy_name_hash,
    };
    write_info(&opts.output_dir, &info)?;

    // Persist the (cleaned) reference sequences (piscem does not retain them) so
    // the selective-alignment step fetches windows matching the indexed k-mers.
    write_refseq_store(&opts.output_dir, std::slice::from_ref(&cleaned), &idx)
        .context("writing reference sequence store")?;

    if !opts.keep_intermediate {
        remove_cdbg_intermediates(&cdbg_prefix);
        if !opts.keep_fixed_fasta {
            let _ = std::fs::remove_file(&cleaned);
        }
    }

    info!("index built: {} references", info.num_refs);
    Ok(info)
}

fn write_info(dir: &Path, info: &IndexInfo) -> Result<()> {
    let path = dir.join(INFO_FILE);
    let json = serde_json::to_string_pretty(info)?;
    std::fs::write(&path, json).with_context(|| format!("writing {}", path.display()))?;
    Ok(())
}

/// First whitespace-delimited token of a FASTA header (matches how cf1/piscem
/// derive reference names).
fn header_name(id: &[u8]) -> &[u8] {
    let end = id
        .iter()
        .position(|b| b.is_ascii_whitespace())
        .unwrap_or(id.len());
    &id[..end]
}

/// The processed reference name: the header up to the first space/tab (salmon's
/// `sepStr`), and — when `gencode` is set — additionally truncated at the first
/// `|` (so GENCODE `ENST...|ENSG...|...` headers reduce to `ENST...`). Used for
/// the cleaned-FASTA headers (→ index names), decoy matching, and the name hash,
/// so all three stay consistent.
fn processed_name(id: &[u8], gencode: bool) -> &[u8] {
    let end = id
        .iter()
        .position(|&b| b == b' ' || b == b'\t')
        .unwrap_or(id.len());
    let name = &id[..end];
    if gencode {
        let pe = name.iter().position(|&b| b == b'|').unwrap_or(name.len());
        &name[..pe]
    } else {
        name
    }
}

/// Read all transcript FASTA(s) and write the reference sequences concatenated
/// in the index's transcript-id order, plus the offset table.
fn write_refseq_store(dir: &Path, transcripts: &[PathBuf], idx: &ReferenceIndex) -> Result<()> {
    use std::collections::HashMap;

    // name -> uppercased sequence bytes
    let mut by_name: HashMap<Vec<u8>, Vec<u8>> = HashMap::new();
    for fasta in transcripts {
        let mut reader = needletail::parse_fastx_file(fasta)
            .with_context(|| format!("opening {}", fasta.display()))?;
        while let Some(rec) = reader.next() {
            let rec = rec.with_context(|| format!("reading {}", fasta.display()))?;
            let name = header_name(rec.id()).to_vec();
            let seq: Vec<u8> = rec.seq().iter().map(|b| b.to_ascii_uppercase()).collect();
            by_name.insert(name, seq);
        }
    }

    let num_refs = idx.num_refs();
    let mut concat: Vec<u8> = Vec::new();
    let mut offsets: Vec<u64> = Vec::with_capacity(num_refs + 1);
    offsets.push(0);
    for tid in 0..num_refs {
        let name = idx.ref_name(tid);
        let seq = by_name.get(name.as_bytes()).ok_or_else(|| {
            anyhow::anyhow!("reference '{name}' from the index was not found in the input FASTA")
        })?;
        concat.extend_from_slice(seq);
        offsets.push(concat.len() as u64);
    }

    std::fs::write(dir.join(REFSEQ_FILE), &concat)
        .with_context(|| format!("writing {REFSEQ_FILE}"))?;
    std::fs::write(dir.join(REFSEQ_OFFSETS_FILE), serde_json::to_vec(&offsets)?)
        .with_context(|| format!("writing {REFSEQ_OFFSETS_FILE}"))?;
    Ok(())
}

fn remove_cdbg_intermediates(cdbg_prefix: &Path) {
    for ext in ["cf_seg", "cf_seq", "json"] {
        let p = cdbg_prefix.with_extension(ext);
        let _ = std::fs::remove_file(p);
    }
}

/// A loaded salmon index: the piscem reference index, the reference sequences,
/// and salmon metadata.
pub struct SalmonIndex {
    inner: ReferenceIndex,
    info: IndexInfo,
    /// concatenated reference sequences (forward strand), transcript-id order
    refseq: Vec<u8>,
    /// offsets into `refseq` (`num_refs + 1` entries)
    ref_offsets: Vec<u64>,
}

impl SalmonIndex {
    /// Load a salmon index from a directory produced by [`build`].
    pub fn load(dir: impl AsRef<Path>) -> Result<Self> {
        let dir = dir.as_ref();
        // Guard: a C++ salmon (pufferfish) index is not readable by salmon 2.0 —
        // the index format changed. Detect it (no sshash `index.ssi`, but
        // pufferfish artifacts present) and give an actionable rebuild message
        // instead of a cryptic parse error from `read_info`.
        if !dir.join(format!("{INDEX_PREFIX}.ssi")).exists()
            && (dir.join("ctable.bin").exists() || dir.join("pre_indexing.log").exists())
        {
            anyhow::bail!(
                "{} looks like a C++ salmon (pufferfish) index, which salmon 2.0 cannot read.\n\
                 The index format changed in 2.0 — rebuild it with `salmon index -t <transcripts> -i <index>`.",
                dir.display()
            );
        }
        let info = read_info(dir)?;
        let index_prefix = dir.join(INDEX_PREFIX);
        // Load the EC table when the index has one (enables pseudoalignment-only mode).
        let inner = ReferenceIndex::load(&index_prefix, info.has_ec_table, false)
            .with_context(|| format!("loading piscem index at {}", index_prefix.display()))?;
        let refseq = std::fs::read(dir.join(REFSEQ_FILE))
            .with_context(|| format!("reading {REFSEQ_FILE}"))?;
        let ref_offsets: Vec<u64> = serde_json::from_slice(
            &std::fs::read(dir.join(REFSEQ_OFFSETS_FILE))
                .with_context(|| format!("reading {REFSEQ_OFFSETS_FILE}"))?,
        )?;
        Ok(Self {
            inner,
            info,
            refseq,
            ref_offsets,
        })
    }

    /// Forward-strand sequence of reference `tid`.
    pub fn ref_seq(&self, tid: u32) -> &[u8] {
        let s = self.ref_offsets[tid as usize] as usize;
        let e = self.ref_offsets[tid as usize + 1] as usize;
        &self.refseq[s..e]
    }

    /// The concatenated forward-strand reference sequences (all transcripts).
    /// Indexed by [`ref_offsets`](Self::ref_offsets); used to build the
    /// reduced-memory GC rank bitvector.
    pub fn refseq_concat(&self) -> &[u8] {
        &self.refseq
    }

    /// Cumulative offsets into [`refseq_concat`](Self::refseq_concat)
    /// (`num_refs + 1` entries; transcript `tid` spans `[off[tid], off[tid+1])`).
    pub fn ref_offsets(&self) -> &[u64] {
        &self.ref_offsets
    }

    /// k-mer length of the index.
    pub fn k(&self) -> usize {
        self.inner.k()
    }

    /// minimizer length of the index.
    pub fn m(&self) -> usize {
        self.inner.m()
    }

    /// number of reference sequences (transcripts).
    pub fn num_refs(&self) -> usize {
        self.inner.num_refs()
    }

    /// name of reference `i`.
    pub fn ref_name(&self, i: usize) -> &str {
        self.inner.ref_name(i)
    }

    /// length of reference `i`.
    pub fn ref_len(&self, i: usize) -> u64 {
        self.inner.ref_len(i)
    }

    /// whether an equivalence-class table is available (pseudoalignment mode).
    pub fn has_ec_table(&self) -> bool {
        self.inner.has_ec_table()
    }

    /// salmon index metadata.
    pub fn info(&self) -> &IndexInfo {
        &self.info
    }

    /// Borrow the underlying piscem index (used by the mapper).
    pub fn inner(&self) -> &ReferenceIndex {
        &self.inner
    }
}

impl salmon_core::RefProvider for SalmonIndex {
    fn num_refs(&self) -> usize {
        self.num_refs()
    }
    fn ref_seq(&self, tid: u32) -> &[u8] {
        self.ref_seq(tid)
    }
    fn is_decoy(&self, tid: u32) -> bool {
        // O(1) range test (matching C++ salmon's `tid >= firstDecoyIndex`). The
        // build guarantees decoys occupy one contiguous block
        // `[first_decoy_index, first_decoy_index + num_decoys)` (it errors out
        // otherwise), so this needs no per-fragment set lookup. The trailing sub-`k`
        // "short" transcripts after the block are NOT decoys (they are 0-count
        // transcripts reported in quant.sf); they carry no k-mers, so `is_decoy` is
        // never consulted for them during mapping regardless.
        let t = tid as usize;
        match self.info.first_decoy_index {
            // Legacy index that recorded `first_decoy_index` but not `num_decoys`
            // (== 0): keep the old behavior — everything from `first_decoy_index`
            // on is a decoy.
            Some(fdi) if self.info.num_decoys == 0 => t >= fdi,
            Some(fdi) => t >= fdi && t < fdi + self.info.num_decoys,
            None => false,
        }
    }
}

fn read_info(dir: &Path) -> Result<IndexInfo> {
    let path = dir.join(INFO_FILE);
    let bytes = std::fs::read(&path).with_context(|| format!("reading {}", path.display()))?;
    let info: IndexInfo =
        serde_json::from_slice(&bytes).with_context(|| format!("parsing {}", path.display()))?;
    Ok(info)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// Write a tiny transcriptome FASTA and return its path.
    fn write_fasta(dir: &Path) -> PathBuf {
        // Three transcripts, each well over k=31, with some shared content so
        // the cDBG has a few branching unitigs.
        let core = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let seqs = [
            ("t1", format!("{core}GGGGCCCCAAAATTTT{core}")),
            ("t2", format!("{core}TTTTAAAACCCCGGGG{core}")),
            (
                "t3",
                format!("GATTACAGATTACAGATTACA{core}CACACACACACACACACACA"),
            ),
        ];
        let path = dir.join("txome.fa");
        let mut f = std::fs::File::create(&path).unwrap();
        for (name, seq) in seqs {
            writeln!(f, ">{name}\n{seq}").unwrap();
        }
        path
    }

    #[test]
    fn build_and_load_roundtrip() {
        let tmp = tempfile::tempdir().unwrap();
        let fasta = write_fasta(tmp.path());
        let out = tmp.path().join("idx");

        let mut opts = IndexBuildOptions::new(vec![fasta], out.clone());
        opts.threads = 1;
        let built = build(&opts).expect("index build");
        assert_eq!(built.k, 31);
        assert_eq!(built.num_refs, 3, "should index 3 transcripts");

        // info.json persisted and the cDBG intermediates were cleaned up.
        assert!(out.join("info.json").exists());
        assert!(!out.join("cdbg.cf_seg").exists());

        let idx = SalmonIndex::load(&out).expect("index load");
        assert_eq!(idx.k(), 31);
        assert_eq!(idx.num_refs(), 3);
        let names: Vec<&str> = (0..idx.num_refs()).map(|i| idx.ref_name(i)).collect();
        assert!(names.contains(&"t1"), "names were {names:?}");
        assert!(idx.ref_len(0) > 31);

        // The reference sequence store round-trips and matches the index lengths.
        for tid in 0..idx.num_refs() as u32 {
            let seq = idx.ref_seq(tid);
            assert_eq!(
                seq.len() as u64,
                idx.ref_len(tid as usize),
                "tid {tid} length"
            );
            assert!(seq.iter().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')));
        }
    }
}
