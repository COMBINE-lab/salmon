//! Salmon index construction and loading.
//!
//! # What an index is, and why one is needed
//!
//! Finding where a 100-base read comes from by comparing it against every
//! transcript would take forever. Instead salmon builds a *search structure* once
//! per transcriptome — the index — and then answers each read's query against it
//! in roughly constant time. Building takes minutes; the payoff is that mapping
//! hundreds of millions of reads becomes feasible.
//!
//! Salmon's index is a piscem (SSHash + contig table) index built over the
//! compacted de Bruijn graph / reference tiling of the transcriptome. In plain
//! terms:
//!
//! * a **k-mer** is a substring of fixed length `k` (31 by default). Every
//!   sequence is a series of overlapping k-mers.
//! * a **de Bruijn graph** links k-mers that overlap by `k-1` bases. Sequence
//!   shared between transcripts becomes *one* path in the graph rather than being
//!   stored repeatedly, which is what makes the structure compact. "Compacted"
//!   means unbranching runs are collapsed into single nodes called *unitigs*.
//! * **SSHash** is a dictionary that answers "which unitig contains this k-mer,
//!   and where?" very fast and in little memory.
//! * the **contig table** maps a unitig back to every transcript/position it came
//!   from, which is how a k-mer hit becomes a list of candidate transcripts.
//!
//! This crate orchestrates the two-stage build that salmon's C++
//! `BuildSalmonIndex` performs, but using the Rust components:
//!
//! 1. [`cf1_rs::cf_build`] constructs the cDBG and writes the
//!    `.cf_seg` / `.cf_seq` / `.json` tiling.
//! 2. [`piscem_rs::index::build::build_index`] builds the SSHash dictionary and
//!    contig table from that tiling.
//!
//! A small `info.json` records the build parameters. [`SalmonIndex`] wraps the
//! loaded [`piscem_rs`] [`ReferenceIndex`] and exposes the accessors the
//! mapper and quantifier need.
//!
//! # Decoys, in one paragraph
//!
//! A decoy is extra sequence (usually the whole genome) indexed alongside the
//! transcripts purely as a sink. Reads from introns, unannotated genes or
//! genomic DNA contamination would otherwise be forced onto whichever transcript
//! they resemble most; with decoys present they align better to the decoy and are
//! discarded. Much of the machinery below exists to guarantee that decoys occupy
//! one *contiguous* block of reference ids, because that turns the
//! per-fragment "is this a decoy?" test into a single range comparison.

use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use cf1_rs::{cf_build, CfInput};
use piscem_rs::index::build::{build_index, BuildConfig};
use piscem_rs::index::reference_index::ReferenceIndex;
use salmon_core::RefSeqs;
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
/// Where `REFSEQ_FILE` is staged, so a rebuild that fails on bad input leaves
/// the existing store untouched.
const REFSEQ_STAGED_FILE: &str = "refseq.bin.partial";
/// Per-reference offsets into `REFSEQ_FILE` (`num_refs + 1` entries).
const REFSEQ_OFFSETS_FILE: &str = "refseq_offsets.json";

/// On-disk salmon index *format* version written into `info.json` by this build.
/// Distinct from `salmon_version` (the software release string): this is bumped
/// only when a change to the index layout or its semantics makes a previously
/// written index unsafe to load with the current code.
///
/// History:
/// - **0** — implicit (no `index_version` field); salmon 2.0.0–2.0.1. The decoy
///   block could be non-contiguous when the decoy set contained a sub-`k`
///   reference (piscem relocated it into the trailing "short" block), which the
///   O(1) `is_decoy` range test would then mis-classify, silently dropping real
///   transcripts. Short/N-untileable decoy handling was also different.
/// - **1** — salmon 2.1.0: sub-`k` / N-untileable decoys are dropped at build so
///   the decoy block is always a single contiguous range (enforced by a
///   build-time contiguity guard), `N` is retained in decoy sequences, and the
///   reference layout `[transcripts][decoys][short transcripts]` is guaranteed.
/// - **2** — sshash-rs 0.7 one-modality minimizer scheme (C++ SSHash v6 port):
///   the `.ssi` dictionary format moved to 5.0 (canonical minimizer at forward
///   density, centre-closest tie-break; header records the build seed and a
///   hasher guard). Format-v1 indices (every salmon release from 2.1.0 up to
///   this format change) cannot be loaded and must be rebuilt. `info.json` now records `sshash_format_version` so
///   future dictionary format breaks are diagnosable from the index dir.
/// - **3** — sshash-rs 0.8: the minimizer MPHF is sux's PHastR (in place of
///   `ph`'s PHast) and sux moved to epserde 0.13, so the `.ssi` dictionary format
///   is 6.0 and the MPHF container 3.0. Format-v2 indices (salmon 2.7.0 – 2.8.0)
///   cannot be loaded and must be rebuilt. Contigs are now numbered by first
///   appearance in the reference tilings, which makes index builds
///   byte-reproducible.
pub const INDEX_FORMAT_VERSION: u32 = 3;

/// Minimum index format version this build can load. Indices below this are
/// rejected with an actionable rebuild message (their decoy / short-transcript
/// layout predates the contiguity guarantee and cannot be loaded safely).
///
/// Refusing to load is deliberate: silently producing subtly wrong abundances is
/// far worse than telling the user to spend ten minutes rebuilding.
pub const MIN_READABLE_INDEX_VERSION: u32 = 3;

/// The sshash `.ssi` dictionary format major version this salmon builds
/// against (recorded in `info.json` as `sshash_format_version`).
pub const SSHASH_FORMAT_VERSION_MAJOR: u32 = 6;

/// Default minimizer length for a given k (used when `m == 0`).
///
/// A *minimizer* is the smallest k-mer-like substring within a window, used as a
/// bucket key so that neighbouring positions share a bucket. Longer minimizers
/// mean fewer, more specific buckets; they must be shorter than `k`.
///
/// Matches piscem's default minimizer length of 19 for the standard k = 31.
/// Falls back to `k/2 + 1` for small k, where 19 would be invalid (m must be < k).
fn default_minimizer_len(k: usize) -> usize {
    if k > 19 {
        19
    } else {
        (k / 2) + 1
    }
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
    ///
    /// Odd is required so a k-mer can never equal its own reverse complement,
    /// which would make the canonical form ambiguous.
    pub k: usize,
    /// minimizer length; `0` selects a default derived from `k`
    pub m: usize,
    /// worker threads; `0` means all available cores
    pub threads: usize,
    /// Ignored since index format v2: the sshash dictionary is always
    /// canonical (a k-mer and its reverse complement share one key — the
    /// salmon/pufferfish convention is now the only modality).
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
    ///
    /// Identical sequences are statistically indistinguishable, so the EM cannot
    /// apportion reads between them; collapsing avoids reporting an arbitrary
    /// split of the same evidence.
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
    ///
    /// Poly-A tails are added post-transcriptionally and are near-identical
    /// across all transcripts, so indexing them creates a large block of shared,
    /// uninformative sequence.
    pub clip_polya: bool,
    /// Directory for sshash's external minimizer-sort scratch. `None` defaults to
    /// a `sshash_tmp` subdirectory of the build-intermediate dir (i.e. `--tmpdir`
    /// when set, else the output dir); `Some(dir)` overrides it, e.g. to place the
    /// sort scratch on a separate/fast disk. Salmon creates the directory before
    /// the build and removes it afterwards unless `keep_intermediate` is set.
    /// (Without this, sshash writes to `sshash_tmp` relative to the current
    /// working directory and leaves the empty directory behind.)
    pub sshash_tmp: Option<PathBuf>,
    /// RAM ceiling, in GiB, for sshash's external minimizer sort (the dominant
    /// build-time memory/disk trade-off). `None` keeps sshash's default (8 GiB);
    /// a smaller value spills to disk sooner (less RAM, more I/O).
    pub ram_limit_gib: Option<usize>,
}

impl IndexBuildOptions {
    /// Construct options with salmon's defaults for the given transcripts and
    /// output directory. Callers then override individual fields.
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
            sshash_tmp: None,
            ram_limit_gib: None,
        }
    }
}

/// Metadata written to `info.json` alongside the piscem index.
///
/// This is what makes an index self-describing: a later quant run reads it to
/// learn the k-mer size, the decoy layout, and which reference FASTA the index
/// was built from (via the hashes).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IndexInfo {
    pub k: usize,
    pub m: usize,
    /// Always `true` since index format v2 (the sshash dictionary is
    /// canonical by construction); kept for readers of older info.json.
    pub canonical: bool,
    pub has_ec_table: bool,
    pub num_refs: usize,
    /// Index of the first decoy reference: references `[first_decoy_index, num_refs)`
    /// are decoys (genome/contamination sequences indexed for selective-alignment
    /// but never quantified). `None` (or absent) means no decoys. salmon requires
    /// all decoy records to appear contiguously at the end of the FASTA.
    ///
    /// `#[serde(default)]` means the field may be missing from an older
    /// `info.json` and will then take its type's default rather than failing to
    /// parse — that is how new fields are added without invalidating old files.
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
    /// Salmon index *format* version (see [`INDEX_FORMAT_VERSION`]). Absent in
    /// indices built by salmon <= 2.0.1, which deserialize to `0` — below
    /// [`MIN_READABLE_INDEX_VERSION`], so those are rejected on load. Distinct
    /// from `salmon_version`, the software release string below.
    #[serde(default)]
    pub index_version: u32,
    /// The sshash `.ssi` dictionary format major version the index was built
    /// with (5 = sshash-rs 0.7 unified scheme). Absent (0) in older indices.
    /// Recorded so a dictionary format break is diagnosable from `info.json`
    /// without poking the binary `.ssi` header.
    #[serde(default)]
    pub sshash_format_version: u32,
    /// whether the index was built with `--keepDuplicates`.
    ///
    /// `Option` rather than a defaulted `bool`: an index built before this was
    /// recorded must read as *unknown*, not as `false`, or a quant run would
    /// report a guess as an observation. `None` propagates into
    /// `meta_info.json`'s `incomplete_meta_info_fields`.
    #[serde(default)]
    pub keep_duplicates: Option<bool>,
    /// salmon version that produced the index
    pub salmon_version: String,
    /// SHA-256/512 of the reference sequences and names, computed to byte-match
    /// salmon's pufferfish `FixFasta` hashes (so downstream provenance tools
    /// recognize a Rust-built index as equivalent for the same FASTA). Empty for
    /// indexes built before these were recorded (`#[serde(default)]`).
    ///
    /// The point is reproducibility: two indexes with the same hashes were built
    /// from the same reference, so results from them are comparable.
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

/// The reference sequence/name SHA-256 and SHA-512 hashes, computed exactly as
/// salmon's pufferfish `FixFasta` does, so the values byte-match a salmon index
/// built from the same FASTA(s):
/// - sequences are hashed in FASTA order, with non-printable bytes removed but
///   **before** uppercasing / non-ACGT replacement / poly-A clipping (i.e. the
///   original-case printable sequence), for every record;
/// - names are the header up to the first space or tab (salmon's `sepStr=" \t"`),
///   hashed for every record whose (printable) sequence is non-empty;
/// - decoy records (named in the decoy list) are routed to dedicated decoy
///   hashers, so a decoy-aware index gets real decoy digests; with no decoys
///   those hashers stay empty ⇒ SHA of "", matching salmon's no-decoy case.
///
/// Every one of those "before"s matters: hashing the *input* rather than the
/// processed sequence is what lets the digest identify the FASTA the user
/// supplied, independently of how salmon later transforms it.
#[derive(Default)]
struct RefDigests {
    seq256: sha2::Sha256,
    name256: sha2::Sha256,
    seq512: sha2::Sha512,
    name512: sha2::Sha512,
    // Decoy sequences/names are routed to their own hashers (salmon's
    // `update_seq_hash(isDecoy,..)` / `update_name_hash`); with no decoys these
    // stay empty → SHA of "" (matches salmon's no-decoy case). Only the SHA-256
    // decoy digests are recorded by salmon, so we only need those.
    decoy_seq256: sha2::Sha256,
    decoy_name256: sha2::Sha256,
    /// the printable bytes of a record that has others
    printable: Vec<u8>,
}

impl RefDigests {
    /// Add one record: its processed name and its sequence as read.
    fn add(&mut self, name: &[u8], seq: &[u8], is_decoy: bool) {
        use sha2::Digest;
        // Keep only C-`isprint` bytes (0x20..=0x7e); preserve original case.
        // This drops line endings and stray control bytes so the digest does
        // not depend on how the FASTA was wrapped.
        let seq = if is_printable(seq) {
            seq
        } else {
            self.printable.clear();
            self.printable
                .extend(seq.iter().copied().filter(|b| (0x20..=0x7e).contains(b)));
            &self.printable
        };
        if is_decoy {
            self.decoy_seq256.update(seq);
        } else {
            self.seq256.update(seq);
            self.seq512.update(seq);
        }
        // Empty-sequence records contribute no name, matching salmon.
        if !seq.is_empty() {
            if is_decoy {
                self.decoy_name256.update(name);
            } else {
                self.name256.update(name);
                self.name512.update(name);
            }
        }
    }

    fn finish(self) -> RefHashes {
        use sha2::Digest;
        // Render a digest as lowercase hex, the conventional textual form.
        let hex = |bytes: &[u8]| -> String {
            let mut s = String::with_capacity(bytes.len() * 2);
            for b in bytes {
                s.push_str(&format!("{b:02x}"));
            }
            s
        };
        RefHashes {
            seq_hash: hex(&self.seq256.finalize()),
            name_hash: hex(&self.name256.finalize()),
            seq_hash512: hex(&self.seq512.finalize()),
            name_hash512: hex(&self.name512.finalize()),
            decoy_seq_hash: hex(&self.decoy_seq256.finalize()),
            decoy_name_hash: hex(&self.decoy_name256.finalize()),
        }
    }
}

/// Whether every byte of `seq` is C-`isprint` (0x20..=0x7e), written to
/// vectorize like [`is_pure_acgt`].
fn is_printable(seq: &[u8]) -> bool {
    seq.chunks(64).all(|block| {
        block
            .iter()
            .fold(0u8, |bad, &b| bad | (b.wrapping_sub(0x20) > 0x5e) as u8)
            == 0
    })
}

/// Hash the reference inputs in a pass of their own, on a thread that overlaps
/// the whole build: SHA-512 costs about twice what the parallel-decoded
/// preprocessing pass does (1.5 Gbases of GENCODE), so it must not hold that
/// pass back, and nothing reads the digests until `info.json` is written.
///
/// This reads every input a second time, which only a regular file allows; see
/// [`HashFeed`] for streams.
fn compute_ref_hashes(
    transcripts: &[PathBuf],
    decoy_names: &ahash::AHashSet<Vec<u8>>,
    gencode: bool,
) -> Result<RefHashes> {
    let mut digests = RefDigests::default();
    for path in transcripts {
        let mut reader = needletail::parse_fastx_file(path)
            .with_context(|| format!("opening {} for hashing", path.display()))?;
        while let Some(rec) = reader.next() {
            let rec = rec.context("reading FASTA record for hashing")?;
            // Header up to the first ' '/'\t' (and first '|' under --gencode).
            let name = processed_name(rec.id(), gencode);
            digests.add(name, &rec.seq(), decoy_names.contains(name));
        }
    }
    Ok(digests.finish())
}

/// The digests of inputs that cannot be read twice (a FIFO, or process
/// substitution such as `-t <(curl ...)`), fed by [`preprocess_fasta`] as it
/// reads them, so that every input is read exactly once.
///
/// It hashes on a thread of its own, from the batches preprocessing hands the
/// writers: a record that leaves preprocessing unchanged is hashed from its
/// batch, and only one that was changed or dropped is sent as read. A stream is
/// decoded serially, which is slower than hashing, so the thread keeps up; the
/// queue lets it fall behind by about 1 GiB at most.
struct HashFeed {
    tx: std::sync::mpsc::SyncSender<(Batch, Vec<HashItem>)>,
}

/// A batch of processed `(name, sequence)` records, shared by the writers.
type Batch = std::sync::Arc<Vec<(Vec<u8>, Vec<u8>)>>;

impl HashFeed {
    /// Start the hashing thread in `scope`; it ends, with the digests, once the
    /// `HashFeed` is dropped.
    fn spawn<'scope>(
        scope: &'scope std::thread::Scope<'scope, '_>,
    ) -> (
        Self,
        std::thread::ScopedJoinHandle<'scope, Result<RefHashes>>,
    ) {
        // Batches hold about 4 Mbases each.
        let (tx, rx) = std::sync::mpsc::sync_channel::<(Batch, Vec<HashItem>)>(256);
        let handle = scope.spawn(move || {
            let mut digests = RefDigests::default();
            for (batch, items) in rx {
                for item in items {
                    match item {
                        HashItem::Kept(i, is_decoy) => {
                            let (name, seq) = &batch[i as usize];
                            digests.add(name, seq, is_decoy);
                        }
                        HashItem::Own(name, seq, is_decoy) => digests.add(&name, &seq, is_decoy),
                    }
                }
            }
            Ok(digests.finish())
        });
        (Self { tx }, handle)
    }

    fn send(&self, batch: Batch, items: Vec<HashItem>) {
        // The thread only stops early by panicking, which the join re-raises.
        let _ = self.tx.send((batch, items));
    }
}

/// The six digests [`RefDigests`] produces, grouped so the function can
/// return them all without a six-element tuple.
struct RefHashes {
    seq_hash: String,
    name_hash: String,
    seq_hash512: String,
    name_hash512: String,
    decoy_seq_hash: String,
    decoy_name_hash: String,
}

/// Outcome of [`preprocess_fasta`]: how many bases were randomized and which
/// transcripts were dropped as exact-sequence duplicates.
struct PreprocessResult {
    /// non-ACGT bases replaced with pseudo-random ACGT
    replaced: u64,
    /// `(retained_name, duplicate_name)` pairs in discovery order. Without
    /// `--keepDuplicates` the duplicate was *not* written to the cleaned FASTA (so
    /// it is absent from the index); with `--keepDuplicates` every transcript is
    /// written and this just records which ones are exact duplicates of each other.
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
    /// the reference sequence store, written alongside the cleaned FASTA
    refseq: StagedRefseq,
}

/// What a [`HashFeed`] is sent for one input record, in input order.
enum HashItem {
    /// The record is entry `.0` of the batch it travels with, unchanged.
    Kept(u32, bool),
    /// The record as read, because it was changed or dropped:
    /// `(name, sequence, is_decoy)`.
    Own(Vec<u8>, Vec<u8>, bool),
}

/// salmon/pufferfish `FixFasta` clips a reference iff its (cleaned) sequence ends
/// in a run of at least this many `A`s; the entire trailing `A` run is then
/// removed (an all-`A` reference is dropped).
///
/// The threshold avoids clipping a transcript that merely happens to end in one
/// or two genuine `A`s.
const POLYA_CLIP_LEN: usize = 10;

/// Length of the longest run of consecutive ACGT (case-insensitive) bases — the
/// longest stretch Cuttlefish can extract k-mers from (it splits the de Bruijn
/// graph on any other base, e.g. `N`). A run of length `L` carries `L - k + 1`
/// k-mers, so a reference is tileable iff its longest ACGT run is `>= k`.
///
/// So a 50 kb genome fragment made of 20-base islands separated by `N` runs is,
/// for indexing purposes, no more useful than a 20-base sequence.
fn longest_acgt_run(seq: &[u8]) -> usize {
    let mut best = 0usize;
    let mut cur = 0usize;
    for &b in seq {
        if matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't') {
            cur += 1;
            best = best.max(cur);
        } else {
            // Any other byte breaks the run.
            cur = 0;
        }
    }
    best
}

/// Read every input FASTA, apply salmon's reference cleaning rules, and write one
/// combined cleaned FASTA for the index build. See [`PreprocessResult`] for what
/// it reports back.
///
/// The rules, applied per record in this order: reject a transcript that follows
/// a decoy, drop exact duplicates, randomize non-ACGT bases (transcripts only),
/// clip poly-A tails, drop untileable decoys.
///
/// Every reference written to the cleaned FASTA is also staged into `refseq`.
/// With a `hash_feed`, every record read, written or not, is also sent to be
/// hashed for `info.json`, so that each input is read exactly once.
#[allow(clippy::too_many_arguments)]
fn preprocess_fasta(
    inputs: &[PathBuf],
    out_path: &Path,
    mut refseq: RefseqWriter,
    threads: usize,
    keep_duplicates: bool,
    decoy_names: &ahash::AHashSet<Vec<u8>>,
    gencode: bool,
    clip_polya: bool,
    k: usize,
    hash_feed: Option<HashFeed>,
) -> Result<PreprocessResult> {
    use std::io::Write as _;
    const B: [u8; 4] = *b"ACGT";
    // Seed of a small deterministic generator (below) used for base replacement.
    // Deterministic on purpose: the same FASTA must always produce the same
    // index, so a real random source would be wrong here.
    let mut x: u64 = 0x9E37_79B9_7F4A_7C15;
    let out_file = std::fs::File::create(out_path)
        .with_context(|| format!("creating {}", out_path.display()))?;

    // The cleaned FASTA and the refseq store are written by threads of their
    // own, fed batches of processed records, so this thread only parses and
    // cleans: the two writes cost as much as everything else here combined.
    const BATCH_BASES: usize = 1 << 22;
    std::thread::scope(|scope| {
        let (fasta_tx, fasta_rx) = std::sync::mpsc::sync_channel::<Batch>(4);
        let (refseq_tx, refseq_rx) = std::sync::mpsc::sync_channel::<Batch>(4);
        let fasta_writer = scope.spawn(move || -> Result<()> {
            let mut out = std::io::BufWriter::with_capacity(1 << 20, out_file);
            for batch in fasta_rx {
                for (name, seq) in batch.iter() {
                    // The record in canonical two-line form (header, then the whole
                    // sequence unwrapped), which is what the downstream builder expects.
                    out.write_all(b">")?;
                    out.write_all(name)?;
                    out.write_all(b"\n")?;
                    out.write_all(seq)?;
                    out.write_all(b"\n")?;
                }
            }
            out.flush()?;
            Ok(())
        });
        let refseq_writer = scope.spawn(move || -> Result<StagedRefseq> {
            for batch in refseq_rx {
                for (name, seq) in batch.iter() {
                    refseq.push(name, seq)?;
                }
            }
            refseq.finish()
        });
        // Hand a batch to both writers and the hasher; `false` once a writer has
        // stopped (failed).
        let send = |batch: Vec<(Vec<u8>, Vec<u8>)>, items: Vec<HashItem>| {
            let batch = std::sync::Arc::new(batch);
            if let Some(feed) = &hash_feed {
                feed.send(batch.clone(), items);
            }
            fasta_tx.send(batch.clone()).is_ok() && refseq_tx.send(batch).is_ok()
        };
        let mut batch = Vec::new();
        let mut hash_items: Vec<HashItem> = Vec::new();
        let feeding = hash_feed.is_some();
        let mut batch_bases = 0usize;
        let mut replaced = 0u64;
        let mut polya_clipped = 0u64;
        let mut polya_dropped: Vec<Vec<u8>> = Vec::new();
        let mut short_decoys_dropped: Vec<Vec<u8>> = Vec::new();
        let mut duplicate_clusters: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
        // Tracks whether a decoy has been seen, to enforce the "decoys last" rule.
        let mut saw_decoy = false;
        // Names from the decoy list not yet found in the input.
        let mut unmatched_decoys = decoy_names.clone();
        // Exact original (pure-ACGT) sequence -> name of the first transcript that
        // carried it. salmon collapses transcripts whose *processed* sequences match;
        // because it randomizes non-ACGT bases with a position-advancing RNG, two
        // identical sequences that contain any non-ACGT base never compare equal and
        // so are never collapsed — we mirror that by only deduplicating pure-ACGT
        // sequences (keyed on the exact original bytes, matching salmon's pre-process
        // XXH64 over the unmodified sequence). Populated in both modes so duplicates
        // are detected for `duplicate_clusters.tsv` even under --keepDuplicates.
        //
        // Keyed by a 128-bit XXH3 digest rather than the sequence itself: holding a
        // copy of every distinct transcript would cost as much memory as the whole
        // transcriptome. salmon itself keys on a 64-bit digest; at 128 bits an
        // accidental collision among even 10^9 sequences has probability ~10^-21.
        let mut seen: ahash::AHashMap<u128, Vec<u8>> = ahash::AHashMap::new();
        let mut names_seen: ahash::AHashSet<Vec<u8>> = ahash::AHashSet::new();
        'records: for path in inputs {
            let mut reader = open_fasta(path, threads)?;
            while let Some(rec) = reader.next() {
                let rec = rec.context("reading FASTA record")?;
                // Processed name: header up to the first space/tab (salmon's sepStr),
                // and the first `|` too under `--gencode`.
                let id = rec.id();
                let name = processed_name(id, gencode);
                let orig = rec.seq();
                // Everything downstream identifies a reference by name: two with one
                // name become rows that gene maps and `quantmerge` cannot tell apart.
                if !names_seen.insert(name.to_vec()) {
                    anyhow::bail!(
                        "reference name {:?} appears more than once in the FASTA input; \
                         names must be unique",
                        String::from_utf8_lossy(name)
                    );
                }
                let is_decoy = decoy_names.contains(name);
                // salmon requires every decoy record to appear, contiguously, at the
                // end of the input — reject a real transcript after any decoy.
                //
                // This is what ultimately guarantees the contiguous decoy id block;
                // failing here is far better than discovering the violation later,
                // when the only symptom is wrong abundances.
                if is_decoy {
                    saw_decoy = true;
                    unmatched_decoys.remove(name);
                } else if saw_decoy {
                    anyhow::bail!(
                        "non-decoy reference {:?} appears after a decoy; decoy records must be \
                     contiguous at the end of the FASTA input",
                        String::from_utf8_lossy(name)
                    );
                }

                // Detect exact-sequence duplicates (pure-ACGT only — salmon hashes the
                // unmodified sequence, and because it randomizes any non-ACGT base
                // per-position, sequences containing one never compare equal). The
                // cluster is recorded for `duplicate_clusters.tsv` in BOTH modes,
                // matching salmon (which reports duplicates even with --keepDuplicates);
                // without --keepDuplicates the duplicate is additionally dropped
                // (collapsed onto the first occurrence).
                let pure_acgt = is_pure_acgt(&orig);
                if pure_acgt {
                    let digest = xxhash_rust::xxh3::xxh3_128(&orig);
                    if let Some(retained_name) = seen.get(&digest) {
                        duplicate_clusters.push((retained_name.clone(), name.to_vec()));
                        if !keep_duplicates {
                            // `continue` skips the write below, so the duplicate never
                            // reaches the index.
                            if feeding {
                                let seq = orig.into_owned();
                                hash_items.push(HashItem::Own(name.to_vec(), seq, is_decoy));
                            }
                            continue;
                        }
                    } else {
                        seen.insert(digest, name.to_vec());
                    }
                }

                let mut seq = orig.into_owned();
                // The sequence as read, kept for the hasher once `seq` changes.
                let mut as_read: Option<Vec<u8>> = None;
                // Transcripts: replace non-ACGT bases with pseudo-random ACGT so the
                // processed sequence is fully k-mer-able (salmon FixFasta behavior).
                // DECOYS: leave non-ACGT bases (e.g. genome N-runs) in place — cf1-rs
                // splits the de Bruijn graph on them natively (and records the N-gaps in
                // the tiling), yielding a far less tangled, smaller, faster-to-build
                // cDBG than seeding spurious k-mers from random replacements across
                // assembly gaps. The refseq store keeps the raw bytes and the aligner
                // encodes N as a mismatch (dna5 code 4), so a read aligning over a decoy
                // N-run simply scores it as a mismatch.
                //
                // A pure-ACGT sequence (nearly every transcript) has nothing to replace.
                if !is_decoy && !pure_acgt {
                    if feeding {
                        as_read = Some(seq.clone());
                    }
                    for b in seq.iter_mut() {
                        if !matches!(*b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't') {
                            // A linear congruential generator: multiply, add, and take
                            // high bits. `wrapping_*` means overflow wraps around
                            // rather than panicking, which is the intended arithmetic.
                            x = x
                                .wrapping_mul(6364136223846793005)
                                .wrapping_add(1442695040888963407);
                            *b = B[((x >> 33) & 3) as usize];
                            replaced += 1;
                        }
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
                    // `rposition` finds the last non-`A`, i.e. where the tail begins.
                    match seq.iter().rposition(|&b| !matches!(b, b'A' | b'a')) {
                        Some(last) => {
                            if feeding && as_read.is_none() {
                                as_read = Some(seq.clone());
                            }
                            seq.truncate(last + 1);
                            polya_clipped += 1;
                        }
                        None => {
                            // All `A`s: drop the reference entirely (not written, not
                            // counted as retained), matching salmon.
                            polya_clipped += 1;
                            polya_dropped.push(name.to_vec());
                            if feeding {
                                let seq = as_read.unwrap_or(seq);
                                hash_items.push(HashItem::Own(name.to_vec(), seq, is_decoy));
                            }
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
                    if feeding {
                        let seq = as_read.unwrap_or(seq);
                        hash_items.push(HashItem::Own(name.to_vec(), seq, is_decoy));
                    }
                    continue;
                }

                if feeding {
                    hash_items.push(match as_read {
                        Some(read) => HashItem::Own(name.to_vec(), read, is_decoy),
                        None => HashItem::Kept(batch.len() as u32, is_decoy),
                    });
                }
                batch_bases += seq.len();
                batch.push((name.to_vec(), seq));
                if batch_bases >= BATCH_BASES {
                    batch_bases = 0;
                    let items = std::mem::take(&mut hash_items);
                    if !send(std::mem::take(&mut batch), items) {
                        // A writer failed; its error is reported below.
                        break 'records;
                    }
                }
            }
        }
        send(batch, hash_items);
        drop((fasta_tx, refseq_tx));
        // A writer panic is re-raised here rather than reported as an error.
        fasta_writer
            .join()
            .unwrap_or_else(|panic| std::panic::resume_unwind(panic))
            .with_context(|| format!("writing {}", out_path.display()))?;
        let refseq = refseq_writer
            .join()
            .unwrap_or_else(|panic| std::panic::resume_unwind(panic))?;
        // Without this, a list that matches nothing builds "successfully" and
        // quantifies the genome as if it were a transcript.
        if !unmatched_decoys.is_empty() {
            let mut names: Vec<_> = unmatched_decoys
                .iter()
                .map(|n| String::from_utf8_lossy(n))
                .collect();
            names.sort();
            anyhow::bail!(
                "{} name(s) in the decoy list match no reference in the FASTA, e.g. {:?}. \
                 Names must be the header up to the first space, with no leading `>`.",
                names.len(),
                &names[..names.len().min(5)]
            );
        }
        Ok(PreprocessResult {
            replaced,
            duplicate_clusters,
            polya_clipped,
            polya_dropped,
            short_decoys_dropped,
            refseq,
        })
    })
}

/// Whether `seq` consists only of `ACGTacgt`.
///
/// Written so the compiler vectorizes it even for baseline x86-64 (a `matches!`
/// on the bytes lowers to a per-byte table test instead): `b & 0xDF` upper-cases
/// exactly the bytes `acgt` and maps no other byte onto `ACGT`, and each 64-byte
/// block is folded branch-free, with one early-exit test per block.
fn is_pure_acgt(seq: &[u8]) -> bool {
    seq.chunks(64).all(|block| {
        block.iter().fold(0u8, |bad, &b| {
            let u = b & 0xDF;
            bad | !((u == b'A') | (u == b'C') | (u == b'G') | (u == b'T')) as u8
        }) == 0
    })
}

/// Parallel gzip decoding of a reference stops paying off here: the single
/// thread parsing the decoded FASTA is then the bottleneck, and more decoders
/// only buffer more decoded data ahead of it.
const MAX_DECODE_THREADS: usize = 4;

/// Open a FASTA (optionally gzipped) for record-by-record reading.
///
/// A gzipped regular file is decoded in parallel (up to `threads` workers, at
/// most [`MAX_DECODE_THREADS`]); anything else goes through needletail's own
/// (serial) decompression.
///
/// Anything but a regular file (a FIFO, or process substitution such as
/// `-t <(curl ...)`) is opened exactly once and handed to needletail, which
/// detects compression from the stream itself: sniffing the magic through a
/// second open would consume the start of the stream.
fn open_fasta(path: &Path, threads: usize) -> Result<Box<dyn needletail::FastxReader>> {
    use std::io::Read as _;
    if !path.metadata().is_ok_and(|m| m.is_file()) {
        let file =
            std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
        return needletail::parse_fastx_reader(file)
            .with_context(|| format!("opening {}", path.display()));
    }
    let mut magic = [0u8; 2];
    let gzipped = std::fs::File::open(path)
        .and_then(|mut f| f.read_exact(&mut magic))
        .is_ok()
        && magic == [0x1f, 0x8b];
    let reader = if gzipped {
        let decoder = rapidgzip_core::Decoder::builder()
            .decoder_threads(threads.clamp(1, MAX_DECODE_THREADS))
            .build()
            .map_err(|e| anyhow::anyhow!("creating gzip decoder: {e}"))?;
        let decoded = decoder
            .open(path)
            .map_err(|e| anyhow::anyhow!("opening {}: {e}", path.display()))?;
        needletail::parse_fastx_reader(decoded)
    } else {
        needletail::parse_fastx_file(path)
    };
    reader.with_context(|| format!("opening {}", path.display()))
}

/// The reference sequence store (`refseq.bin`) as staged during preprocessing,
/// in the order the index is expected to number the references.
struct StagedRefseq {
    /// reference names in expected index order
    names: Vec<Vec<u8>>,
    /// `names.len() + 1` offsets into `refseq.bin`
    offsets: Vec<u64>,
}

/// Writes `refseq.bin` while the references are being preprocessed, so the
/// build never has to read the cleaned FASTA back to produce it.
///
/// The index numbers references in cleaned-FASTA order, except that references
/// shorter than `k` (which carry no k-mers and so never enter the de Bruijn
/// graph) are moved, in order, to the end. The writer mirrors that: long
/// references stream straight to disk, short ones are held back and appended
/// on [`finish`](Self::finish). Whether the built index really used this order
/// is checked against its reference names afterwards ([`finish_refseq_store`]).
struct RefseqWriter {
    file: std::io::BufWriter<std::fs::File>,
    k: usize,
    names: Vec<Vec<u8>>,
    offsets: Vec<u64>,
    /// sub-`k` references, held back to be appended last
    shorts: Vec<(Vec<u8>, Vec<u8>)>,
    /// upper-casing scratch, reused across records
    upper: Vec<u8>,
}

impl RefseqWriter {
    fn create(path: &Path, k: usize) -> Result<Self> {
        let file =
            std::fs::File::create(path).with_context(|| format!("creating {}", path.display()))?;
        Ok(Self {
            file: std::io::BufWriter::with_capacity(1 << 20, file),
            k,
            names: Vec::new(),
            offsets: vec![0],
            shorts: Vec::new(),
            upper: Vec::new(),
        })
    }

    /// Append one processed reference (stored upper-cased, as the aligner reads
    /// it: soft-masking is irrelevant to alignment).
    fn push(&mut self, name: &[u8], seq: &[u8]) -> Result<()> {
        if seq.len() < self.k {
            self.shorts.push((name.to_vec(), seq.to_ascii_uppercase()));
            return Ok(());
        }
        self.upper.clear();
        self.upper.extend_from_slice(seq);
        self.upper.make_ascii_uppercase();
        self.write(name.to_vec(), seq.len())
    }

    fn write(&mut self, name: Vec<u8>, len: usize) -> Result<()> {
        use std::io::Write as _;
        self.file
            .write_all(&self.upper)
            .with_context(|| format!("writing {REFSEQ_FILE}"))?;
        let end = self.offsets.last().copied().unwrap_or(0) + len as u64;
        self.offsets.push(end);
        self.names.push(name);
        Ok(())
    }

    fn finish(mut self) -> Result<StagedRefseq> {
        use std::io::Write as _;
        for (name, seq) in std::mem::take(&mut self.shorts) {
            self.upper = seq;
            let len = self.upper.len();
            self.write(name, len)?;
        }
        self.file
            .flush()
            .with_context(|| format!("writing {REFSEQ_FILE}"))?;
        Ok(StagedRefseq {
            names: self.names,
            offsets: self.offsets,
        })
    }
}

/// Complete the reference sequence store staged during preprocessing: if the
/// built index numbers its references exactly as staged (the normal case), only
/// the offset table remains to be written; otherwise the store is rebuilt in the
/// index's order from the cleaned FASTA.
fn finish_refseq_store(
    dir: &Path,
    staged: StagedRefseq,
    cleaned: &Path,
    idx: &ReferenceIndex,
) -> Result<()> {
    let in_order = staged.names.len() == idx.num_refs()
        && staged
            .names
            .iter()
            .enumerate()
            .all(|(tid, name)| idx.ref_name(tid).as_bytes() == name.as_slice());
    if !in_order {
        warn!(
            "index reference order differs from the preprocessed order; \
             rebuilding {REFSEQ_FILE} from the cleaned FASTA"
        );
        let _ = std::fs::remove_file(dir.join(REFSEQ_STAGED_FILE));
        return write_refseq_store(dir, std::slice::from_ref(&cleaned.to_path_buf()), idx);
    }
    std::fs::rename(dir.join(REFSEQ_STAGED_FILE), dir.join(REFSEQ_FILE))
        .with_context(|| format!("writing {REFSEQ_FILE}"))?;
    std::fs::write(
        dir.join(REFSEQ_OFFSETS_FILE),
        serde_json::to_vec(&staged.offsets)?,
    )
    .with_context(|| format!("writing {REFSEQ_OFFSETS_FILE}"))?;
    Ok(())
}

/// Read a decoy-names file (one name per line; blank lines ignored) into a set,
/// matching how `name` is extracted in [`preprocess_fasta`] (no `>` prefix).
///
/// A set (rather than a list) because membership is tested once per FASTA record
/// and the decoy list can hold thousands of contig names.
fn read_decoy_names(path: &Path) -> Result<ahash::AHashSet<Vec<u8>>> {
    let bytes =
        std::fs::read(path).with_context(|| format!("reading decoys file {}", path.display()))?;
    let mut set = ahash::AHashSet::new();
    for line in bytes.split(|&b| b == b'\n') {
        // trim trailing CR and surrounding whitespace
        //
        // Done by hand on raw bytes (rather than with a string method) because
        // reference names are not guaranteed to be valid UTF-8. The `[rest @ ..,
        // last]` syntax is slice pattern matching: it binds the last element and
        // everything before it.
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
///
/// Unconditional output matters for pipelines: a missing file is an error to
/// handle, an empty one is not.
fn write_duplicate_clusters(dir: &Path, clusters: &[(Vec<u8>, Vec<u8>)]) -> Result<()> {
    use std::io::Write as _;
    let path = dir.join("duplicate_clusters.tsv");
    let mut f = std::io::BufWriter::new(
        std::fs::File::create(&path).with_context(|| format!("creating {}", path.display()))?,
    );
    f.write_all(b"RetainedRef\tDuplicateRef\n")?;
    for (retained, dup) in clusters {
        // Written as raw bytes rather than formatted text, since names need not
        // be valid UTF-8.
        f.write_all(retained)?;
        f.write_all(b"\t")?;
        f.write_all(dup)?;
        f.write_all(b"\n")?;
    }
    f.flush()?;
    Ok(())
}

/// Build a salmon index from a transcriptome FASTA into `opts.output_dir`.
///
/// The whole build in order: validate options, read the decoy list, preprocess
/// the FASTA, build the cDBG (stage 1), build the SSHash index (stage 2), locate
/// and verify the decoy block, write the reference-sequence store and `info.json`,
/// and clean up intermediates.
pub fn build(opts: &IndexBuildOptions) -> Result<IndexInfo> {
    anyhow::ensure!(!opts.transcripts.is_empty(), "no transcript files provided");
    // tiny-dict panics above 31, after the build has written files; odd is
    // required so no k-mer equals its own reverse complement.
    anyhow::ensure!(
        opts.k >= 1 && opts.k <= 31 && opts.k % 2 == 1,
        "k must be odd in [1, 31], got {}",
        opts.k
    );
    anyhow::ensure!(
        opts.m < opts.k,
        "the minimizer length must be less than k, got m={}, k={}",
        opts.m,
        opts.k
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

    // The reference digests cover the *input* FASTA and nothing reads them
    // until `info.json` is written, so they are computed on a thread of their
    // own that overlaps the whole build. Regular files are hashed in a pass of
    // their own, as before; when any input is a stream, which can be read only
    // once, preprocessing feeds the hashing instead. The scope joins the thread
    // even if the build fails.
    let streamed = opts
        .transcripts
        .iter()
        .any(|p| !p.metadata().is_ok_and(|m| m.is_file()));
    std::thread::scope(|scope| {
        let (hash_feed, ref_hashes) = if streamed {
            let (feed, handle) = HashFeed::spawn(scope);
            (Some(feed), handle)
        } else {
            let handle =
                scope.spawn(|| compute_ref_hashes(&opts.transcripts, &decoy_names, opts.gencode));
            (None, handle)
        };

        // Preprocess the references into the cleaned FASTA that feeds both the cDBG
        // build and the refseq store (the hashes below use the original input):
        // transcripts get non-ACGT bases replaced with pseudo-random ACGT (salmon
        // FixFasta behavior); decoys keep their non-ACGT bases (e.g. genome N-runs) so
        // Cuttlefish splits the graph on them instead of tiling random replacements.
        let cleaned = intermediate_dir.join("cleaned_refs.fa");
        let refseq_staged = opts.output_dir.join(REFSEQ_STAGED_FILE);
        let refseq = RefseqWriter::create(&refseq_staged, opts.k)?;
        let pre = preprocess_fasta(
            &opts.transcripts,
            &cleaned,
            refseq,
            threads,
            opts.keep_duplicates,
            &decoy_names,
            opts.gencode,
            opts.clip_polya,
            opts.k,
            hash_feed,
        )
        .context("preprocessing reference FASTA (non-ACGT replacement, dedup)")
        .inspect_err(|_| {
            let _ = std::fs::remove_file(&refseq_staged);
        })?;
        // The input is good, so the rebuild is committed: from here the directory
        // holds a mix of old and new files until the final `info.json` is written.
        // An index with no `info.json` refuses to load, which is the safe failure.
        match std::fs::remove_file(opts.output_dir.join(INFO_FILE)) {
            Err(e) if e.kind() != std::io::ErrorKind::NotFound => {
                return Err(e).context("removing the previous info.json");
            }
            _ => {}
        }
        // Everything below reports what the preprocessing did. Silent transformation
        // of a user's reference would be indefensible, so each rule that fired is
        // logged, with the destructive ones (dropped references) as warnings.
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
        if !pre.duplicate_clusters.is_empty() {
            if opts.keep_duplicates {
                info!(
                    "detected {} exact-sequence-duplicate transcripts; retained them \
                     (--keepDuplicates) and recorded them in duplicate_clusters.tsv",
                    pre.duplicate_clusters.len()
                );
            } else {
                info!(
                    "removed {} transcripts that were exact sequence duplicates \
                     (use --keepDuplicates to retain them)",
                    pre.duplicate_clusters.len()
                );
            }
        }
        // Always emit duplicate_clusters.tsv (header + one row per duplicate), matching
        // salmon — downstream tooling expects it, and with --keepDuplicates it is the
        // record of which retained transcripts are exact duplicates of each other.
        write_duplicate_clusters(&opts.output_dir, &pre.duplicate_clusters)
            .context("writing duplicate_clusters.tsv")?;

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

        // sshash sorts its minimizer tuples through an on-disk scratch directory.
        // Left to its default it writes to `sshash_tmp` relative to the *current
        // working directory* and leaves the (emptied) directory behind. Anchor it
        // under our build-intermediate dir instead (a `sshash_tmp` subdir of
        // `--tmpdir`, or the output dir when unset), or wherever the caller asked,
        // so the scratch is co-located with the other intermediates and we can
        // remove it cleanly afterwards.
        let sshash_tmp = opts
            .sshash_tmp
            .clone()
            .unwrap_or_else(|| intermediate_dir.join("sshash_tmp"));
        std::fs::create_dir_all(&sshash_tmp)
            .with_context(|| format!("creating sshash scratch dir {}", sshash_tmp.display()))?;

        // Stage 2: piscem SSHash + contig-table index over the tiling.
        info!("building piscem index (m={})", m);
        let config = BuildConfig {
            input_prefix: cdbg_prefix.clone(),
            output_prefix: index_prefix.clone(),
            k: opts.k,
            m,
            build_ec_table: opts.build_ec_table,
            num_threads: opts.threads,
            // Fixed seed: the hash function must be reproducible, or two builds of
            // the same reference would differ.
            seed: 1,
            mphf_build: Default::default(),
            emit_tiny: None,
            tmp_dir: Some(sshash_tmp.clone()),
            ram_limit_gib: opts.ram_limit_gib,
        };
        build_index(&config).context("piscem index construction failed")?;

        // Load once to capture reference count and confirm the index is usable.
        // Better to fail here, during the build, than at the start of a quant run.
        let idx = ReferenceIndex::load(&index_prefix, opts.build_ec_table, false)
            .context("loading freshly built index")?;
        // Reference seq/name hashes (salmon-FixFasta-compatible), over the input FASTA.
        let h = ref_hashes
            .join()
            .map_err(|_| anyhow::anyhow!("reference hashing thread panicked"))?
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
            // Which reference ids ended up being decoys, in index order.
            let decoy_tids: Vec<usize> = (0..n)
                .filter(|&t| decoy_names.contains(idx.ref_name(t).as_bytes()))
                .collect();
            match decoy_tids.first().copied() {
                // Named decoys existed but none survived preprocessing.
                None => (None, 0usize),
                Some(first) => {
                    let count = decoy_tids.len();
                    // contiguous iff the decoy tids are exactly first..first+count
                    let contiguous = decoy_tids.iter().enumerate().all(|(i, &t)| t == first + i);
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
            // Everything after the decoy block is a short (sub-k) transcript.
            let n_short = idx.num_refs() - (fdi + num_decoys);
            info!(
                "{num_decoys} decoy reference(s) (first decoy index {fdi}); \
                 {n_short} short (sub-k, 0-count) transcript(s) recorded after the decoys"
            );
        }
        let info = IndexInfo {
            k: opts.k,
            m,
            canonical: true,
            has_ec_table: idx.has_ec_table(),
            num_refs: idx.num_refs(),
            first_decoy_index,
            num_decoys,
            index_version: INDEX_FORMAT_VERSION,
            sshash_format_version: SSHASH_FORMAT_VERSION_MAJOR,
            keep_duplicates: Some(opts.keep_duplicates),
            // `env!` reads an environment variable at compile time; Cargo sets this
            // one from Cargo.toml, so the recorded version is always this build's.
            salmon_version: env!("CARGO_PKG_VERSION").to_string(),
            seq_hash: h.seq_hash,
            name_hash: h.name_hash,
            seq_hash512: h.seq_hash512,
            name_hash512: h.name_hash512,
            decoy_seq_hash: h.decoy_seq_hash,
            decoy_name_hash: h.decoy_name_hash,
        };

        // Persist the (cleaned) reference sequences (piscem does not retain them) so
        // the selective-alignment step fetches windows matching the indexed k-mers.
        finish_refseq_store(&opts.output_dir, pre.refseq, &cleaned, &idx)
            .context("writing reference sequence store")?;
        // Last, because it is what marks the index complete.
        write_info(&opts.output_dir, &info)?;

        if !opts.keep_intermediate {
            remove_cdbg_intermediates(&cdbg_prefix);
            if !opts.keep_fixed_fasta {
                // `let _ =` deliberately ignores the result: failing to delete a
                // temporary file must not fail an otherwise successful build.
                let _ = std::fs::remove_file(&cleaned);
            }
            // sshash removes its own scratch files but leaves the directory. Prune
            // it: `remove_dir` only succeeds on an empty directory, so a caller-
            // supplied override that is shared / pre-populated is left untouched.
            let _ = std::fs::remove_dir(&sshash_tmp);
        }

        info!("index built: {} references", info.num_refs);
        Ok(info)
    })
}

/// Serialize [`IndexInfo`] to `info.json` (pretty-printed, since users read it).
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
///
/// That consistency is the point: if decoy matching used a different name rule
/// than the index, a decoy would silently fail to be recognized.
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
///
/// **Why store the sequences at all.** piscem's index holds k-mers, not
/// sequence, but selective alignment has to compare a read against the actual
/// reference bases, and the bias models need the local sequence context. Storing
/// one concatenated blob plus an offset table is the cheapest structure that
/// gives O(1) access to any transcript: `refseq[offsets[tid]..offsets[tid+1]]`.
///
/// Order matters: the sequences are written in the *index's* id order, not the
/// FASTA's, because that is how every later lookup indexes them.
fn write_refseq_store(dir: &Path, transcripts: &[PathBuf], idx: &ReferenceIndex) -> Result<()> {
    use std::collections::HashMap;

    // name -> uppercased sequence bytes
    //
    // Uppercased because FASTA files use lower case for soft-masked (repetitive)
    // regions, which is irrelevant to alignment and would otherwise have to be
    // handled at every comparison.
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
    // `num_refs + 1` entries: the extra final offset is the end of the last
    // sequence, so every transcript's span is `offsets[tid]..offsets[tid+1]`
    // with no special case for the last one.
    let mut offsets: Vec<u64> = Vec::with_capacity(num_refs + 1);
    offsets.push(0);
    for tid in 0..num_refs {
        let name = idx.ref_name(tid);
        let seq = by_name.get(name.as_bytes()).ok_or_else(|| {
            anyhow::anyhow!("reference '{name}' from the index was not found in the input FASTA")
        })?;
        anyhow::ensure!(
            seq.len() as u64 == idx.ref_len(tid),
            "reference '{name}' has length {} in the FASTA but {} in the index \
             (duplicate reference names?)",
            seq.len(),
            idx.ref_len(tid)
        );
        concat.extend_from_slice(seq);
        offsets.push(concat.len() as u64);
    }

    std::fs::write(dir.join(REFSEQ_FILE), &concat)
        .with_context(|| format!("writing {REFSEQ_FILE}"))?;
    std::fs::write(dir.join(REFSEQ_OFFSETS_FILE), serde_json::to_vec(&offsets)?)
        .with_context(|| format!("writing {REFSEQ_OFFSETS_FILE}"))?;
    Ok(())
}

/// Delete the cDBG tiling files, which are only needed as input to stage 2.
/// Failures are ignored: leftover temporaries are untidy, not incorrect.
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
    /// concatenated reference sequences (forward strand), transcript-id order.
    /// Empty when the index was loaded with `load_refseq == false` (sketch mode
    /// without sequence-dependent bias correction): the raw nucleotides are only
    /// needed for selective-alignment extension and seq/GC bias, so skipping the
    /// read avoids holding a multi-GB decoy genome resident.
    refseq: Vec<u8>,
    /// the same sequences at two bits per base, in place of `refseq` under
    /// [`RefseqLoad::Packed`]
    packed: PackedRefseq,
    /// offsets into `refseq` (`num_refs + 1` entries); empty iff `!refseq_loaded`
    ref_offsets: Vec<u64>,
    /// whether the reference sequence bytes were loaded
    refseq_loaded: bool,
}

/// How [`SalmonIndex::load_with_opts`] holds the reference sequences.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RefseqLoad {
    /// Not at all: nothing will read them.
    None,
    /// One byte per base, for callers of [`SalmonIndex::ref_seq`] (bias models).
    Bytes,
    /// Two bits per base, for mapping, which reads them through
    /// [`RefProvider::ref_view`](salmon_core::RefProvider::ref_view).
    Packed,
}

/// `refseq.bin` at two bits per base; see [`salmon_core::pack_bases`].
#[derive(Default)]
struct PackedRefseq {
    bits: Vec<u8>,
    other: Vec<salmon_core::OtherRun>,
    /// For each reference, the range of `other` that overlaps it (a run that
    /// crosses from one reference into the next is in both). Looking the runs
    /// up by binary search on every `ref_view` cost 5% of a genome-decoy run,
    /// which calls it once per seed.
    runs_of: Vec<(u32, u32)>,
}

impl PackedRefseq {
    /// Pack the file directly into its final buffer, so that one byte per base
    /// is never in memory; block by block across threads, so loading takes no
    /// longer than reading the file whole.
    fn load(path: &Path) -> Result<Self> {
        use std::io::{Read, Seek, SeekFrom};
        /// bases per block; a multiple of 4 so that blocks start on a packed byte
        const BLOCK: usize = 1 << 22;
        let len = std::fs::metadata(path)
            .with_context(|| format!("reading {REFSEQ_FILE}"))?
            .len() as usize;
        let mut bits = vec![0u8; len.div_ceil(4)];
        let threads = std::thread::available_parallelism().map_or(1, |n| n.get().min(8));
        let per_thread = (bits.len().div_ceil(threads)).max(1);
        let parts: Vec<Vec<salmon_core::OtherRun>> = std::thread::scope(|scope| {
            let workers: Vec<_> = bits
                .chunks_mut(per_thread)
                .enumerate()
                .map(|(t, out)| {
                    scope.spawn(move || -> Result<Vec<salmon_core::OtherRun>> {
                        let first = (t * per_thread * 4) as u64;
                        let mut f = std::fs::File::open(path)?;
                        f.seek(SeekFrom::Start(first))?;
                        let total = (len - first as usize).min(out.len() * 4);
                        let mut buf = vec![0u8; BLOCK.min(total)];
                        let mut runs = Vec::new();
                        let mut done = 0;
                        for o in out.chunks_mut(BLOCK / 4) {
                            let n = (total - done).min(BLOCK);
                            f.read_exact(&mut buf[..n])?;
                            salmon_core::pack_bases(&buf[..n], o, first + done as u64, &mut runs);
                            done += n;
                        }
                        Ok(runs)
                    })
                })
                .collect();
            workers
                .into_iter()
                .map(|w| w.join().expect("refseq loader panicked"))
                .collect::<Result<Vec<_>>>()
        })
        .with_context(|| format!("reading {REFSEQ_FILE}"))?;
        // A run can cross from one block or thread into the next, which only
        // leaves it in two pieces: lookups take the first that covers a position.
        let other = parts.concat();
        Ok(Self {
            bits,
            other,
            runs_of: Vec::new(),
        })
    }
}

/// Load just the per-reference forward sequences from an index directory, in
/// transcript-id order (decoys included), without loading the (multi-GB) SSHash
/// dictionary. Reads only `refseq.bin` + `refseq_offsets.json`. Used to feed
/// bias models a requant pass (e.g. `--deterministic`) the reference sequences
/// the index already stores, so the user need not re-supply a transcript FASTA.
pub fn load_ref_seqs(dir: impl AsRef<Path>) -> Result<RefSeqs> {
    let dir = dir.as_ref();
    let refseq =
        std::fs::read(dir.join(REFSEQ_FILE)).with_context(|| format!("reading {REFSEQ_FILE}"))?;
    let offsets: Vec<u64> = serde_json::from_slice(
        &std::fs::read(dir.join(REFSEQ_OFFSETS_FILE))
            .with_context(|| format!("reading {REFSEQ_OFFSETS_FILE}"))?,
    )?;
    // The on-disk layout is already one concatenated blob plus an offset table,
    // which is exactly what `RefSeqs` holds, so both parts move in as they are:
    // no copy of the bases, and no per-transcript allocation. `from_concatenated`
    // validates the table against the blob, turning a corrupt index into an error
    // here rather than a panicking slice at the first bias lookup.
    RefSeqs::from_concatenated(refseq, offsets)
        .with_context(|| format!("reference sequences in {}", dir.display()))
}

impl SalmonIndex {
    /// Load a salmon index from a directory produced by [`build`], including the
    /// reference sequence bytes (needed for selective alignment and seq/GC bias).
    pub fn load(dir: impl AsRef<Path>) -> Result<Self> {
        Self::load_with_opts(dir, RefseqLoad::Bytes)
    }

    /// Load a salmon index, optionally skipping the reference sequence bytes.
    ///
    /// Pass `load_refseq == false` when neither selective-alignment extension nor
    /// seq/GC bias correction will run (e.g. `--sketch` without `--seqBias`/
    /// `--gcBias`): [`ref_seq`](Self::ref_seq) / [`refseq_concat`](Self::refseq_concat)
    /// then must not be called, but the index loads without the (potentially
    /// multi-GB, with a genome decoy) reference sequence resident.
    pub fn load_with_opts(dir: impl AsRef<Path>, load_refseq: RefseqLoad) -> Result<Self> {
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
        // Reject an index whose format predates the decoy / short-transcript
        // contiguity guarantee (salmon <= 2.0.1, which wrote no `index_version`
        // and so deserializes to 0). Such an index can mis-classify decoys and
        // silently drop transcripts, so force a rebuild rather than load it.
        if info.index_version < MIN_READABLE_INDEX_VERSION {
            // Old indices may not record a version string either; say something
            // useful in both cases.
            let built_by = if info.salmon_version.is_empty() {
                "salmon <= 2.0.1".to_string()
            } else {
                format!("salmon {}", info.salmon_version)
            };
            anyhow::bail!(
                "{} was built by {} (index format v{}), which this salmon cannot read \
                 (requires index format v{}+).\n\
                 Index format v2 used an older sshash dictionary format, v1 predates \
                 the unified sshash minimizer scheme, and formats before v1 could \
                 mis-handle decoy and short (sub-k) transcript references; rebuild it \
                 with `salmon index -t <transcripts> [-d <decoys>] -i {}`.",
                dir.display(),
                built_by,
                info.index_version,
                MIN_READABLE_INDEX_VERSION,
                dir.display(),
            );
        }
        // The mirror gate: an index from a FUTURE format would otherwise slip
        // past the >= MIN check and die deep inside the sshash loader with a
        // raw I/O error. Say what actually happened instead.
        if info.index_version > INDEX_FORMAT_VERSION {
            anyhow::bail!(
                "{} was built by salmon {} (index format v{}), which is newer than this \
                 salmon understands (format v{}).\n\
                 Upgrade salmon, or rebuild the index with this version: \
                 `salmon index -t <transcripts> [-d <decoys>] -i {}`.",
                dir.display(),
                info.salmon_version,
                info.index_version,
                INDEX_FORMAT_VERSION,
                dir.display(),
            );
        }

        let index_prefix = dir.join(INDEX_PREFIX);
        // Load the EC table when the index has one (enables pseudoalignment-only mode).
        let inner = ReferenceIndex::load(&index_prefix, info.has_ec_table, false)
            .with_context(|| format!("loading piscem index at {}", index_prefix.display()))?;
        let refseq = match load_refseq {
            RefseqLoad::Bytes => std::fs::read(dir.join(REFSEQ_FILE))
                .with_context(|| format!("reading {REFSEQ_FILE}"))?,
            _ => Vec::new(),
        };
        let packed = match load_refseq {
            RefseqLoad::Packed => PackedRefseq::load(&dir.join(REFSEQ_FILE))?,
            _ => PackedRefseq::default(),
        };
        let ref_offsets: Vec<u64> = if load_refseq == RefseqLoad::None {
            Vec::new()
        } else {
            serde_json::from_slice(
                &std::fs::read(dir.join(REFSEQ_OFFSETS_FILE))
                    .with_context(|| format!("reading {REFSEQ_OFFSETS_FILE}"))?,
            )?
        };
        // A store that disagrees with the index would panic at the first slice or,
        // worse, align reads against the wrong bases.
        if load_refseq != RefseqLoad::None {
            let store_len = std::fs::metadata(dir.join(REFSEQ_FILE))
                .with_context(|| format!("reading {REFSEQ_FILE}"))?
                .len();
            anyhow::ensure!(
                ref_offsets.len() == inner.num_refs() + 1 && ref_offsets.last() == Some(&store_len),
                "{} is damaged: {REFSEQ_FILE} does not match the index. Rebuild it with \
                 `salmon index`.",
                dir.display()
            );
        }
        let mut packed = packed;
        if !packed.bits.is_empty() {
            let other = &packed.other;
            packed.runs_of = ref_offsets
                .windows(2)
                .map(|w| {
                    let lo = other.partition_point(|r| r.start + u64::from(r.len) <= w[0]);
                    let hi = other.partition_point(|r| r.start < w[1]);
                    (lo as u32, hi.max(lo) as u32)
                })
                .collect();
        }
        Ok(Self {
            inner,
            info,
            refseq,
            packed,
            ref_offsets,
            refseq_loaded: load_refseq != RefseqLoad::None,
        })
    }

    /// Whether the reference sequence bytes were loaded (see [`load_with_opts`](Self::load_with_opts)).
    pub fn refseq_loaded(&self) -> bool {
        self.refseq_loaded
    }

    /// Forward-strand sequence of reference `tid`.
    ///
    /// Panics unless the index was loaded with [`RefseqLoad::Bytes`]; callers
    /// in sketch-without-bias mode must not reach the alignment/bias paths that
    /// call this, and mapping reads through `ref_view`.
    pub fn ref_seq(&self, tid: u32) -> &[u8] {
        debug_assert!(
            self.refseq_loaded,
            "ref_seq called on an index loaded without reference sequence"
        );
        // Slice the concatenated blob using the offset table.
        let s = self.ref_offsets[tid as usize] as usize;
        let e = self.ref_offsets[tid as usize + 1] as usize;
        &self.refseq[s..e]
    }

    /// The concatenated forward-strand reference sequences (all transcripts).
    /// Indexed by [`ref_offsets`](Self::ref_offsets); used to build the
    /// reduced-memory GC rank bitvector. Empty if loaded without the reference
    /// sequence (`load_refseq == false`).
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

// Implementing the shared trait is what lets `salmon-map` align against this
// index without depending on this crate.
impl salmon_core::RefProvider for SalmonIndex {
    fn num_refs(&self) -> usize {
        self.num_refs()
    }
    #[inline]
    fn ref_view(&self, tid: u32) -> salmon_core::RefView<'_> {
        if self.packed.bits.is_empty() {
            return salmon_core::RefView::Bytes(self.ref_seq(tid));
        }
        let start = self.ref_offsets[tid as usize];
        let end = self.ref_offsets[tid as usize + 1];
        let (lo, hi) = self.packed.runs_of[tid as usize];
        salmon_core::RefView::Packed {
            bits: &self.packed.bits,
            start,
            len: (end - start) as usize,
            other: &self.packed.other[lo as usize..hi as usize],
        }
    }
    fn is_decoy(&self, tid: u32) -> bool {
        // O(1) range test (matching C++ salmon's `tid >= firstDecoyIndex`). The
        // build guarantees decoys occupy one contiguous block
        // `[first_decoy_index, first_decoy_index + num_decoys)` (it errors out
        // otherwise), so this needs no per-fragment set lookup. The trailing sub-`k`
        // "short" transcripts after the block are NOT decoys (they are 0-count
        // transcripts reported in quant.sf); they carry no k-mers, so `is_decoy` is
        // never consulted for them during mapping regardless.
        //
        // This function runs once per candidate placement — billions of times in a
        // large run — which is why the whole contiguity apparatus above exists.
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

/// Read and parse `info.json` from an index directory.
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

    /// Two records named `t1` gave two `t1` rows in `quant.sf`, which gene maps
    /// and `quantmerge` cannot tell apart. Refused whether or not the sequences
    /// match, and however the build treats duplicate sequences.
    #[test]
    fn a_repeated_transcript_name_is_refused() {
        let random = |mut x: u64| -> String {
            (0..300)
                .map(|_| {
                    x = x
                        .wrapping_mul(6364136223846793005)
                        .wrapping_add(1442695040888963407);
                    ['A', 'C', 'G', 'T'][(x >> 33) as usize & 3]
                })
                .collect()
        };
        let (a, b) = (random(1), random(2));
        for (keep_duplicates, second) in [(false, &b), (false, &a), (true, &a)] {
            let tmp = tempfile::tempdir().unwrap();
            let fasta = tmp.path().join("t.fa");
            std::fs::write(&fasta, format!(">t1\n{a}\n>t2\n{b}\n>t1 again\n{second}\n")).unwrap();
            let mut opts = IndexBuildOptions::new(vec![fasta], tmp.path().join("idx"));
            opts.threads = 1;
            opts.keep_duplicates = keep_duplicates;
            let err = build(&opts).map(|_| ()).unwrap_err();
            let msg = format!("{err:#}");
            assert!(msg.contains("\"t1\""), "{msg}");
        }
    }

    /// The end-to-end contract: a build produces a loadable index whose names,
    /// lengths and stored sequences all agree, and which leaves no intermediates
    /// behind.
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
        // sshash's scratch defaults under the output dir and is pruned after the
        // build (no stray `sshash_tmp` left behind — the bug #1045 addressed).
        assert!(
            !out.join("sshash_tmp").exists(),
            "sshash scratch dir should be pruned after a default build"
        );

        let idx = SalmonIndex::load(&out).expect("index load");
        assert_eq!(idx.k(), 31);
        assert_eq!(idx.num_refs(), 3);
        let names: Vec<&str> = (0..idx.num_refs()).map(|i| idx.ref_name(i)).collect();
        assert!(names.contains(&"t1"), "names were {names:?}");
        assert!(idx.ref_len(0) > 31);

        // The reference sequence store round-trips and matches the index lengths.
        // A mismatch here would mean every alignment reads the wrong bases.
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

    /// A rebuild that fails on bad input must leave the existing index working.
    #[test]
    fn failed_rebuild_keeps_the_existing_index() {
        let tmp = tempfile::tempdir().unwrap();
        let fasta = write_fasta(tmp.path());
        let out = tmp.path().join("idx");
        let mut opts = IndexBuildOptions::new(vec![fasta], out.clone());
        opts.threads = 1;
        build(&opts).expect("index build");
        let before = std::fs::read(out.join(REFSEQ_FILE)).unwrap();

        opts.transcripts = vec![tmp.path().join("typo.fa")];
        assert!(build(&opts).is_err());

        assert_eq!(std::fs::read(out.join(REFSEQ_FILE)).unwrap(), before);
        assert!(!out.join(REFSEQ_STAGED_FILE).exists());
        assert_eq!(SalmonIndex::load(&out).expect("index load").num_refs(), 3);
    }

    /// A store cut short (by anything) is refused at load rather than panicking
    /// at the first read that lands in the missing bases.
    #[test]
    fn a_truncated_reference_store_is_refused() {
        let tmp = tempfile::tempdir().unwrap();
        let fasta = write_fasta(tmp.path());
        let out = tmp.path().join("idx");
        let mut opts = IndexBuildOptions::new(vec![fasta], out.clone());
        opts.threads = 1;
        build(&opts).expect("index build");
        std::fs::write(out.join(REFSEQ_FILE), b"").unwrap();
        let err = SalmonIndex::load(&out)
            .err()
            .expect("must refuse")
            .to_string();
        assert!(err.contains("damaged"), "{err}");
    }

    /// The sshash scratch dir honors an explicit override and the
    /// `keep_intermediate` flag: with `keep_intermediate` it survives at the
    /// overridden location; without it, it is pruned. In neither case does a
    /// `sshash_tmp` appear relative to the current working directory.
    #[test]
    fn sshash_scratch_override_and_keep() {
        let tmp = tempfile::tempdir().unwrap();
        let fasta = write_fasta(tmp.path());
        let scratch = tmp.path().join("fastdisk").join("scratch");

        // keep_intermediate => the overridden scratch dir is retained.
        let mut opts = IndexBuildOptions::new(vec![fasta.clone()], tmp.path().join("idx_keep"));
        opts.threads = 1;
        opts.sshash_tmp = Some(scratch.clone());
        opts.keep_intermediate = true;
        build(&opts).expect("index build (keep)");
        assert!(
            scratch.is_dir(),
            "overridden scratch dir should be kept with keep_intermediate"
        );

        // default (no keep) => the same overridden dir is pruned afterwards.
        let scratch2 = tmp.path().join("fastdisk2").join("scratch");
        let mut opts = IndexBuildOptions::new(vec![fasta], tmp.path().join("idx_prune"));
        opts.threads = 1;
        opts.sshash_tmp = Some(scratch2.clone());
        build(&opts).expect("index build (prune)");
        assert!(
            !scratch2.exists(),
            "overridden scratch dir should be pruned without keep_intermediate"
        );
    }

    /// Build ONE index that simultaneously exercises every short/decoy edge case
    /// and verifies they compose correctly:
    ///   * a sub-`k` **transcript** — KEPT (piscem relocates it to the trailing
    ///     untiled block, reported at 0 reads); the regression that motivated the
    ///     decoy-contiguity guard,
    ///   * a sub-`k` **decoy** — DROPPED at preprocess (no `k`-mer to seed),
    ///   * an N-fragmented **decoy** longer than `k` whose every ACGT run is `< k`
    ///     — DROPPED (longest ACGT run governs tileability, not raw length),
    ///   * a **decoy** with an internal N-run flanked by `>= k` ACGT runs — KEPT,
    ///     with its `N` bytes preserved verbatim in the reference store.
    ///
    /// Then asserts the surviving decoy is one contiguous block so `is_decoy`
    /// stays an O(1) range test.
    ///
    /// One index rather than four, because the interaction is the risk: each rule
    /// is easy in isolation and the failure mode was them shifting each other's
    /// reference ids.
    #[test]
    fn build_composes_short_transcript_short_and_n_decoys() {
        use salmon_core::RefProvider; // brings `is_decoy` into scope
                                      // Lengths below are relative to the default k = 31.
                                      // Distinct ACGT content with no consecutive repeats, so nothing collapses
                                      // as a duplicate and no trailing poly-A run triggers clipping.
        let acgt = |len: usize, seed: usize| -> String {
            const B: [u8; 4] = *b"ACGT";
            (0..len).map(|i| B[(i * 3 + seed) % 4] as char).collect()
        };

        let t_long = acgt(80, 0); // normal transcript: KEPT, tiled
        let t_short = acgt(20, 1); // sub-k transcript: KEPT (0 count), relocated to tail
                                   // decoy with N flanked by two >= k ACGT runs: KEPT, N preserved verbatim
        let d_long = format!("{}NNNNN{}", acgt(40, 2), acgt(40, 3));
        let d_short = acgt(25, 4); // sub-k decoy: DROPPED
                                   // len 46 (> k) but every ACGT run (20, 25) is < k: DROPPED
        let d_nfrag = format!("{}N{}", acgt(20, 5), acgt(25, 6));

        let tmp = tempfile::tempdir().unwrap();
        let fasta = tmp.path().join("gentrome.fa");
        {
            let mut f = std::fs::File::create(&fasta).unwrap();
            // Transcripts first; decoys contiguous at the end (salmon requirement).
            for (n, s) in [
                ("t_long", &t_long),
                ("t_short", &t_short),
                ("d_long", &d_long),
                ("d_short", &d_short),
                ("d_nfrag", &d_nfrag),
            ] {
                writeln!(f, ">{n}\n{s}").unwrap();
            }
        }
        let decoys = tmp.path().join("decoys.txt");
        std::fs::write(&decoys, "d_long\nd_short\nd_nfrag\n").unwrap();

        let out = tmp.path().join("idx");
        let mut opts = IndexBuildOptions::new(vec![fasta], out.clone());
        opts.threads = 1;
        opts.decoys = Some(decoys);
        let built = build(&opts).expect("index build");

        // Dropped d_short + d_nfrag; kept t_long, t_short, d_long.
        assert_eq!(
            built.num_refs, 3,
            "kept t_long + t_short + d_long; dropped d_short + d_nfrag"
        );
        assert_eq!(
            built.first_decoy_index,
            Some(1),
            "the one surviving decoy forms a contiguous block at index 1"
        );
        assert_eq!(built.num_decoys, 1);

        let idx = SalmonIndex::load(&out).expect("index load");
        assert_eq!(idx.num_refs(), 3);
        let names: Vec<&str> = (0..idx.num_refs()).map(|i| idx.ref_name(i)).collect();
        assert!(names.contains(&"t_long"), "names {names:?}");
        assert!(
            names.contains(&"t_short"),
            "sub-k transcript must be kept; names {names:?}"
        );
        assert!(
            names.contains(&"d_long"),
            "tileable N-decoy must be kept; names {names:?}"
        );
        assert!(
            !names.contains(&"d_short"),
            "sub-k decoy must be dropped; names {names:?}"
        );
        assert!(
            !names.contains(&"d_nfrag"),
            "N-fragmented decoy must be dropped; names {names:?}"
        );

        // Decoy classification is exactly {d_long}, via the O(1) range test.
        for (i, n) in names.iter().enumerate() {
            assert_eq!(
                idx.is_decoy(i as u32),
                *n == "d_long",
                "is_decoy({n}) misclassified"
            );
        }

        // The sub-k transcript is reported at its real (untiled) length.
        let t_short_tid = names.iter().position(|n| *n == "t_short").unwrap();
        assert_eq!(idx.ref_len(t_short_tid), t_short.len() as u64);

        // The surviving decoy keeps its N bytes verbatim (decoys are NOT
        // N-replaced, unlike transcripts).
        let d_long_tid = names.iter().position(|n| *n == "d_long").unwrap() as u32;
        let seq = idx.ref_seq(d_long_tid);
        assert_eq!(seq.len(), d_long.len(), "decoy stored at full length");
        assert_eq!(
            seq.iter().filter(|&&b| b == b'N').count(),
            5,
            "all 5 decoy Ns preserved in the reference store"
        );
    }

    /// A decoy list that names something absent from the FASTA is a mismatch to
    /// report, not a reason to index the genome as transcripts.
    #[test]
    fn unmatched_decoy_names_are_refused() {
        let tmp = tempfile::tempdir().unwrap();
        let fasta = write_fasta(tmp.path());
        let decoys = tmp.path().join("decoys.txt");
        std::fs::write(&decoys, ">t3\n").unwrap();
        let mut opts = IndexBuildOptions::new(vec![fasta], tmp.path().join("idx"));
        opts.threads = 1;
        opts.decoys = Some(decoys);
        let err = format!("{:#}", build(&opts).unwrap_err());
        assert!(err.contains(">t3"), "{err}");
    }

    /// An index too old to trust must be refused with a message that says what to
    /// do, rather than loaded and quietly mis-quantified.
    #[test]
    fn rejects_pre_2_1_index() {
        let tmp = tempfile::tempdir().unwrap();
        let fasta = write_fasta(tmp.path());
        let out = tmp.path().join("idx");
        let mut opts = IndexBuildOptions::new(vec![fasta], out.clone());
        opts.threads = 1;
        build(&opts).expect("index build");

        // A freshly built index loads and records the current format version.
        let idx = SalmonIndex::load(&out).expect("current index loads");
        assert_eq!(idx.info().index_version, INDEX_FORMAT_VERSION);

        // Simulate a salmon <= 2.0.1 index: strip the `index_version` field (it
        // deserializes to 0 via `#[serde(default)]`) and set an old software
        // version string. Loading it must fail with an actionable message.
        let info_path = out.join(INFO_FILE);
        let mut info: serde_json::Value =
            serde_json::from_slice(&std::fs::read(&info_path).unwrap()).unwrap();
        info.as_object_mut().unwrap().remove("index_version");
        info["salmon_version"] = serde_json::Value::String("2.0.1".into());
        std::fs::write(&info_path, serde_json::to_vec_pretty(&info).unwrap()).unwrap();

        let err = match SalmonIndex::load(&out) {
            Ok(_) => panic!("pre-2.1 index must be rejected"),
            Err(e) => e,
        };
        // Assert on the *content* of the message: an error the user cannot act on
        // is barely better than a wrong answer.
        let msg = format!("{err:#}");
        assert!(msg.contains("salmon 2.0.1"), "message was: {msg}");
        assert!(msg.contains("rebuild"), "message was: {msg}");
    }

    #[test]
    fn rejects_future_index() {
        // An index claiming a NEWER format than this build understands must be
        // rejected up front with an "upgrade salmon" message, not by a raw
        // sshash I/O error from deep inside the loader.
        let tmp = tempfile::tempdir().unwrap();
        let fasta = write_fasta(tmp.path());
        let out = tmp.path().join("idx");
        let mut opts = IndexBuildOptions::new(vec![fasta], out.clone());
        opts.threads = 1;
        build(&opts).expect("index build");

        let info_path = out.join(INFO_FILE);
        let mut info: serde_json::Value =
            serde_json::from_slice(&std::fs::read(&info_path).unwrap()).unwrap();
        info["index_version"] = serde_json::Value::from(INDEX_FORMAT_VERSION + 1);
        info["salmon_version"] = serde_json::Value::String("99.0.0".into());
        std::fs::write(&info_path, serde_json::to_vec_pretty(&info).unwrap()).unwrap();

        let err = match SalmonIndex::load(&out) {
            Ok(_) => panic!("future-format index must be rejected"),
            Err(e) => e,
        };
        let msg = format!("{err:#}");
        assert!(msg.contains("newer"), "message was: {msg}");
        assert!(msg.contains("Upgrade salmon"), "message was: {msg}");
    }

    #[test]
    fn rejects_v1_index() {
        // A format-v1 index (salmon releases before this format change) predates the unified sshash
        // minimizer scheme: its .ssi is unreadable by sshash-lib 0.7, so the
        // info.json gate must reject it up front with a rebuild message.
        let tmp = tempfile::tempdir().unwrap();
        let fasta = write_fasta(tmp.path());
        let out = tmp.path().join("idx");
        let mut opts = IndexBuildOptions::new(vec![fasta], out.clone());
        opts.threads = 1;
        build(&opts).expect("index build");

        let info_path = out.join(INFO_FILE);
        let mut info: serde_json::Value =
            serde_json::from_slice(&std::fs::read(&info_path).unwrap()).unwrap();
        info["index_version"] = serde_json::Value::from(1u32);
        info["salmon_version"] = serde_json::Value::String("2.4.1".into());
        std::fs::write(&info_path, serde_json::to_vec_pretty(&info).unwrap()).unwrap();

        let err = match SalmonIndex::load(&out) {
            Ok(_) => panic!("format-v1 index must be rejected"),
            Err(e) => e,
        };
        let msg = format!("{err:#}");
        assert!(msg.contains("salmon 2.4.1"), "message was: {msg}");
        assert!(msg.contains("rebuild"), "message was: {msg}");
        assert!(msg.contains("index format v1"), "message was: {msg}");
    }
}
