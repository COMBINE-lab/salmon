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
use tracing::info;

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
fn compute_ref_hashes(transcripts: &[PathBuf]) -> Result<RefHashes> {
    use sha2::{Digest, Sha256, Sha512};
    let mut seq256 = Sha256::new();
    let mut name256 = Sha256::new();
    let mut seq512 = Sha512::new();
    let mut name512 = Sha512::new();
    // No decoy support yet: these stay empty → SHA of "" (matches salmon).
    let decoy256 = Sha256::new();
    let decoy_name256 = Sha256::new();

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
            seq256.update(&seq);
            seq512.update(&seq);
            if !seq.is_empty() {
                // Header up to the first ' ' or '\t'.
                let id = rec.id();
                let end = id
                    .iter()
                    .position(|&b| b == b' ' || b == b'\t')
                    .unwrap_or(id.len());
                let name = &id[..end];
                name256.update(name);
                name512.update(name);
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
        decoy_seq_hash: hex(&decoy256.finalize()),
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
fn preprocess_fasta(inputs: &[PathBuf], out_path: &Path) -> Result<u64> {
    use std::io::Write as _;
    const B: [u8; 4] = *b"ACGT";
    let mut x: u64 = 0x9E37_79B9_7F4A_7C15;
    let mut out = std::io::BufWriter::new(
        std::fs::File::create(out_path).with_context(|| format!("creating {}", out_path.display()))?,
    );
    let mut replaced = 0u64;
    for path in inputs {
        let mut reader = needletail::parse_fastx_file(path)
            .with_context(|| format!("opening {}", path.display()))?;
        while let Some(rec) = reader.next() {
            let rec = rec.context("reading FASTA record")?;
            out.write_all(b">")?;
            out.write_all(rec.id())?;
            out.write_all(b"\n")?;
            let mut seq = rec.seq().into_owned();
            for b in seq.iter_mut() {
                if !matches!(*b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't') {
                    x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                    *b = B[((x >> 33) & 3) as usize];
                    replaced += 1;
                }
            }
            out.write_all(&seq)?;
            out.write_all(b"\n")?;
        }
    }
    out.flush()?;
    Ok(replaced)
}

/// Build a salmon index from a transcriptome FASTA into `opts.output_dir`.
pub fn build(opts: &IndexBuildOptions) -> Result<IndexInfo> {
    anyhow::ensure!(!opts.transcripts.is_empty(), "no transcript files provided");
    anyhow::ensure!(opts.k >= 1 && opts.k <= 63 && opts.k % 2 == 1, "k must be odd in [1, 63]");

    std::fs::create_dir_all(&opts.output_dir)
        .with_context(|| format!("creating index dir {}", opts.output_dir.display()))?;

    let threads = effective_threads(opts.threads);
    let m = if opts.m == 0 { default_minimizer_len(opts.k) } else { opts.m };
    let cdbg_prefix = opts.output_dir.join(CDBG_PREFIX);
    let index_prefix = opts.output_dir.join(INDEX_PREFIX);

    // Replace non-ACGT bases up front (cf1-rs/packed-seq can't index N); the
    // cleaned FASTA feeds the cDBG build and the refseq store, while the hashes
    // below use the original input.
    let cleaned = opts.output_dir.join("cleaned_refs.fa");
    let n_replaced = preprocess_fasta(&opts.transcripts, &cleaned)
        .context("preprocessing reference FASTA (non-ACGT replacement)")?;
    if n_replaced > 0 {
        info!("replaced {n_replaced} non-ACGT base(s) with random bases");
    }

    // Stage 1: compacted de Bruijn graph / reference tiling.
    info!("building compacted dBG (k={}, threads={})", opts.k, threads);
    let cf = cf_build()
        .input(CfInput::Files(vec![cleaned.clone()]))
        .output_prefix(cdbg_prefix.clone())
        .k(opts.k)
        .threads(threads)
        .call()
        .context("cf1-rs cDBG construction failed")?;
    info!(
        "cDBG: {} unitigs, {} distinct k-mers",
        cf.unitig_count, cf.vertex_count
    );

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
    let h = compute_ref_hashes(&opts.transcripts).context("hashing reference sequences/names")?;
    let info = IndexInfo {
        k: opts.k,
        m,
        canonical: opts.canonical,
        has_ec_table: idx.has_ec_table(),
        num_refs: idx.num_refs(),
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
        let _ = std::fs::remove_file(&cleaned);
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
    std::fs::write(
        dir.join(REFSEQ_OFFSETS_FILE),
        serde_json::to_vec(&offsets)?,
    )
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
    // TODO: decoy tracking (info.json) — for now no references are decoys.
}

fn read_info(dir: &Path) -> Result<IndexInfo> {
    let path = dir.join(INFO_FILE);
    let bytes = std::fs::read(&path).with_context(|| format!("reading {}", path.display()))?;
    let info: IndexInfo = serde_json::from_slice(&bytes)
        .with_context(|| format!("parsing {}", path.display()))?;
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
            ("t3", format!("GATTACAGATTACAGATTACA{core}CACACACACACACACACACA")),
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
            assert_eq!(seq.len() as u64, idx.ref_len(tid as usize), "tid {tid} length");
            assert!(seq.iter().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')));
        }
    }
}
