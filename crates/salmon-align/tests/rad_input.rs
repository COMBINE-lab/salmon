//! RAD-input quantification (`quantify_rad`): write a small salmon RAD with a
//! known mapping structure, quantify it, and check counts, mass conservation,
//! profile handling, and thread-count invariance.

use std::path::Path;

use salmon_align::{quantify_rad, AlignQuantOptions, FldPolicy};
use salmon_rad::{
    frag_map_type, FragmentChunkBuf, RadHit, RadOutputWriter, RadProfile, SalmonBulkContext,
    FRAG_LEN_UNPAIRED,
};

/// Write a salmon RAD file. Each fragment is a list of `(tid, frag_len)` hits;
/// `frag_len == None` marks an orphan/single-end placement. If `baked_log_pmf` is
/// `Some`, it is baked into the header (so the reader uses it instead of deriving
/// the FLD from unique fragments).
fn write_rad(
    path: &Path,
    ref_names: &[&str],
    ref_lengths: &[u32],
    profile: RadProfile,
    fragments: &[Vec<(u32, Option<u16>, i32)>],
    baked_log_pmf: Option<&[f64]>,
) {
    write_rad_codec(
        path,
        ref_names,
        ref_lengths,
        profile,
        fragments,
        baked_log_pmf,
        salmon_rad::ChunkCodec::None,
    );
}

#[allow(clippy::too_many_arguments)]
/// Provenance for RADs these tests write, standing in for a mapping pass.
fn prov() -> salmon_rad::WriterProvenance {
    salmon_rad::WriterProvenance {
        mapping_type: salmon_rad::MappingType::Mapping,
        index: None,
        source_programs: vec![],
    }
}

fn write_rad_codec(
    path: &Path,
    ref_names: &[&str],
    ref_lengths: &[u32],
    profile: RadProfile,
    fragments: &[Vec<(u32, Option<u16>, i32)>],
    baked_log_pmf: Option<&[f64]>,
    codec: salmon_rad::ChunkCodec,
) {
    let mut w = RadOutputWriter::create(
        path,
        ref_names,
        ref_lengths,
        true,
        profile,
        1024,
        codec,
        &prov(),
    )
    .unwrap();
    if let Some(pmf) = baked_log_pmf {
        w.set_frag_length_dist(pmf);
    }
    let ctx: SalmonBulkContext = *w.context();
    let mut cb = FragmentChunkBuf::with_capacity_codec(4096, codec);
    for hits in fragments.iter() {
        let proper = hits.iter().all(|(_, fl, _)| fl.is_some());
        let frag_type = if proper {
            frag_map_type::MAPPED_PAIR
        } else {
            frag_map_type::SINGLE
        };
        let rad_hits = hits
            .iter()
            // inward (FR) proper pairs: read1 forward, mate reverse — compatible
            // with the `IU` library type the tests quantify under.
            .map(|&(tid, fl, score)| RadHit {
                tid,
                is_fw: true,
                mate_fw: false,
                pos: 0,
                frag_len: fl.unwrap_or(FRAG_LEN_UNPAIRED),
                score,
            })
            .collect();
        let rec = salmon_rad::SalmonBulkRecord::new(frag_type, rad_hits);
        cb.write(&rec, &ctx).unwrap();
        // emit multiple chunks (real output flushes one chunk per minibatch)
        if cb.nrec() >= 64 {
            w.append_chunk_bytes(&cb.take_bytes().unwrap()).unwrap();
        }
    }
    if cb.nrec() > 0 {
        w.append_chunk_bytes(&cb.take_bytes().unwrap()).unwrap();
    }
    w.finalize().unwrap();
}

fn opts_for(out: &Path) -> AlignQuantOptions {
    let mut o = AlignQuantOptions::new(out.join("unused.rad"), out.to_path_buf());
    o.lib_type = "IU".to_string();
    o.fld_mean = 200.0;
    o.fld_sd = 20.0;
    o
}

/// A RAD whose baked FLD (mean ~500) is deliberately far from the length its
/// records imply (200), so which distribution a run used is unambiguous.
/// 200 fragments on `t0`, 300 on `t1`, all unique proper pairs.
fn write_baked_fld_rad(rad: &Path) {
    let names = ["t0", "t1"];
    let lens = [4000u32, 4000];
    let mut frags: Vec<Vec<(u32, Option<u16>, i32)>> = Vec::new();
    for (tid, n) in [(0u32, 200), (1, 300)] {
        for _ in 0..n {
            frags.push(vec![(tid, Some(200), 0)]);
        }
    }
    // Unnormalized Gaussian centered at 500; `from_log_pmf` normalizes. The
    // reserved slot length is 1024 (see `write_rad`).
    let baked: Vec<f64> = (0..1024)
        .map(|l| {
            let d = (l as f64 - 500.0) / 20.0;
            -0.5 * d * d
        })
        .collect();
    write_rad(rad, &names, &lens, RadProfile::Sketch, &frags, Some(&baked));
}

#[test]
/// With no ambiguity there is nothing to infer, so the counts must come back
/// exactly as given.
fn unique_mappings_recover_exact_counts() {
    let tmp = tempfile::tempdir().unwrap();
    let rad = tmp.path().join("maps.rad");
    let names = ["t0", "t1", "t2"];
    let lens = [1000u32, 1000, 1000];

    // 100 frags -> t0, 200 -> t1, 300 -> t2, all unique proper pairs (frag_len 200).
    let mut frags: Vec<Vec<(u32, Option<u16>, i32)>> = Vec::new();
    for (tid, n) in [(0u32, 100), (1, 200), (2, 300)] {
        for _ in 0..n {
            frags.push(vec![(tid, Some(200), 0)]);
        }
    }
    write_rad(&rad, &names, &lens, RadProfile::Sketch, &frags, None);

    let out = tmp.path().join("quant");
    std::fs::create_dir_all(&out).unwrap();
    let res = quantify_rad(&opts_for(&out), &rad).expect("quantify_rad");

    assert_eq!(res.names, vec!["t0", "t1", "t2"]);
    assert_eq!(res.num_processed, 600);
    assert_eq!(res.num_mapped, 600);
    // unique mappings -> exact per-transcript counts
    for (i, expect) in [100.0, 200.0, 300.0].into_iter().enumerate() {
        assert!(
            (res.counts[i] - expect).abs() < 1.0,
            "t{i}: got {} expected {expect}",
            res.counts[i]
        );
    }
    // mass conservation
    let sum: f64 = res.counts.iter().sum();
    assert!((sum - 600.0).abs() < 1.0, "Σcounts = {sum}");
    // quant.sf written
    assert!(out.join("quant.sf").exists());
}

#[test]
/// Ambiguous fragments must follow the unique evidence, which is the core EM
/// behaviour reaching through the RAD path.
fn multimapping_resolves_toward_higher_abundance() {
    let tmp = tempfile::tempdir().unwrap();
    let rad = tmp.path().join("maps.rad");
    let names = ["t0", "t1"];
    let lens = [1000u32, 1000];

    let mut frags: Vec<Vec<(u32, Option<u16>, i32)>> = Vec::new();
    for _ in 0..300 {
        frags.push(vec![(0, Some(200), 0)]); // unique to t0
    }
    for _ in 0..100 {
        frags.push(vec![(1, Some(200), 0)]); // unique to t1
    }
    for _ in 0..100 {
        frags.push(vec![(0, Some(200), 0), (1, Some(200), 0)]); // ambiguous t0/t1
    }
    write_rad(&rad, &names, &lens, RadProfile::Sketch, &frags, None);

    let out = tmp.path().join("quant");
    std::fs::create_dir_all(&out).unwrap();
    let res = quantify_rad(&opts_for(&out), &rad).expect("quantify_rad");

    assert_eq!(res.num_mapped, 500);
    let sum: f64 = res.counts.iter().sum();
    assert!((sum - 500.0).abs() < 1.0, "Σcounts = {sum}");
    // EM should split the 100 ambiguous fragments ~3:1 toward t0 (its higher
    // unique evidence), so t0 ends well above its 300 unique floor.
    assert!(res.counts[0] > res.counts[1], "{:?}", res.counts);
    assert!(
        res.counts[0] > 360.0 && res.counts[0] < 390.0,
        "t0 = {}",
        res.counts[0]
    );
}

/// Chunk compression is transparent + lossless: quantifying an lz4- or
/// zstd-compressed RAD must yield exactly the same counts as the uncompressed
/// RAD (the reader decompresses chunks in its reader thread). This exercises the
/// full `quantify_rad` -> `ParallelRadReader` decompression path.
#[test]
/// Chunk compression is a storage choice and must not change a single number.
fn compressed_rad_quantifies_identically() {
    use salmon_rad::ChunkCodec;
    let names = ["t0", "t1", "t2"];
    let lens = [1000u32, 1000, 1000];
    let mut frags: Vec<Vec<(u32, Option<u16>, i32)>> = Vec::new();
    for (tid, n) in [(0u32, 120), (1, 250), (2, 70)] {
        for _ in 0..n {
            frags.push(vec![(tid, Some(200), 0)]);
        }
    }
    // a few ambiguous fragments so the EM (not just unique counts) is exercised
    for _ in 0..60 {
        frags.push(vec![(0, Some(200), 0), (1, Some(200), 0)]);
    }

    let quant_with = |codec: ChunkCodec| -> Vec<f64> {
        let tmp = tempfile::tempdir().unwrap();
        let rad = tmp.path().join("maps.rad");
        write_rad_codec(&rad, &names, &lens, RadProfile::Sketch, &frags, None, codec);
        let out = tmp.path().join("q");
        std::fs::create_dir_all(&out).unwrap();
        quantify_rad(&opts_for(&out), &rad)
            .expect("quantify_rad")
            .counts
    };

    let base = quant_with(ChunkCodec::None);
    for codec in [ChunkCodec::Lz4, ChunkCodec::Zstd] {
        let got = quant_with(codec);
        assert_eq!(
            got, base,
            "codec {codec:?} must give identical counts to uncompressed"
        );
    }
}

/// Hand-write a piscem `map-bulk`-style RAD (read tag `frag_map_type`, alignment
/// triple `compressed_ori_ref`/`pos`/`frag_len`, no `frag_name_hash`) and confirm
/// salmon detects the piscem-bulk profile and quantifies it.
#[test]
/// A RAD written by piscem rather than salmon must quantify too — the format is
/// shared, and this is what makes the two tools interoperate.
fn piscem_bulk_input_quantifies() {
    use libradicl::header::{RadHeader, RadPrelude};
    use libradicl::rad_types::{RadIntId, RadType, TagDesc, TagSection, TagSectionLabel, TagValue};
    use libradicl::writers::RadFileWriter;
    use std::io::Write as _;

    let tmp = tempfile::tempdir().unwrap();
    let rad = tmp.path().join("piscem.rad");
    let names = vec!["t0".to_string(), "t1".to_string()];
    let lens = vec![1000u32, 1000];

    // prelude: piscem bulk tag signature (no frag_name_hash).
    let (file_tags, file_tag_map) = TagSection::from_tag_values(
        TagSectionLabel::FileTags,
        &[("ref_lengths", TagValue::ArrayU32(lens.clone()))],
    );
    let mut read_tags = TagSection::new_with_label(TagSectionLabel::ReadTags);
    read_tags.add_tag_desc(TagDesc {
        name: "frag_map_type".into(),
        typeid: RadType::Int(RadIntId::U8),
    });
    let mut aln_tags = TagSection::new_with_label(TagSectionLabel::AlignmentTags);
    for n in ["compressed_ori_ref", "pos", "frag_len"] {
        aln_tags.add_tag_desc(TagDesc {
            name: n.into(),
            typeid: RadType::Int(if n == "frag_len" {
                RadIntId::U16
            } else {
                RadIntId::U32
            }),
        });
    }
    let prelude = RadPrelude {
        hdr: RadHeader {
            is_paired: 1,
            ref_count: names.len() as u64,
            ref_names: names.clone(),
            num_chunks: 0,
        },
        file_tags,
        read_tags,
        aln_tags,
    };
    let mut fw = RadFileWriter::new(
        std::fs::File::create(&rad).unwrap(),
        &prelude,
        &file_tag_map,
    )
    .unwrap();

    // 300 frags -> t0, 100 -> t1 (unique, inward FR proper pairs), in one chunk.
    let mut fragments: Vec<u32> = vec![0u32; 300];
    fragments.extend(std::iter::repeat_n(1u32, 100));
    let mut body: Vec<u8> = Vec::new();
    for &tid in &fragments {
        body.extend(1u32.to_le_bytes()); // na = 1
        body.push(frag_map_type::MAPPED_PAIR); // frag_map_type
                                               // compressed_ori_ref: tid | is_fw bit (FR inward; mate_fw=0)
        let ori = tid | 0x8000_0000;
        body.extend(ori.to_le_bytes());
        body.extend(0u32.to_le_bytes()); // pos
        body.extend(200u16.to_le_bytes()); // frag_len
    }
    let mut chunk = Vec::new();
    chunk.extend(((body.len() + 8) as u32).to_le_bytes()); // nbytes (incl header)
    chunk.extend((fragments.len() as u32).to_le_bytes()); // nrec
    chunk.extend_from_slice(&body);
    fw.write_chunk_bytes(&chunk).unwrap();
    let mut f = fw.finalize().unwrap();
    f.flush().unwrap();

    let out = tmp.path().join("quant");
    std::fs::create_dir_all(&out).unwrap();
    let res = quantify_rad(&opts_for(&out), &rad).expect("quantify_rad (piscem)");

    assert_eq!(res.names, names);
    assert_eq!(res.num_processed, 400);
    assert_eq!(res.num_mapped, 400);
    assert!(
        (res.counts[0] - 300.0).abs() < 1.0,
        "t0 = {}",
        res.counts[0]
    );
    assert!(
        (res.counts[1] - 100.0).abs() < 1.0,
        "t1 = {}",
        res.counts[1]
    );
}

/// Regression test for the worker-drain race: a `quantify_rad` run must never
/// silently return zero (or partial) counts.
///
/// `thread_count_invariant` cannot catch this. It asserts only that the three
/// runs *agree*, so a fault that hits every run looks like a pass — with the
/// race window forced open, all three returned `processed=0` and the test went
/// green. This test asserts the absolute total instead, so a lost chunk fails
/// regardless of how many runs it affects.
///
/// The race: libradicl's producer pushes every meta-chunk and only then sets
/// `done`. A worker whose `pop()` came up empty just before that push would
/// observe `done` and break, abandoning the queue.
#[test]
/// Every input fragment must be accounted for; a silent loss would understate
/// abundances with no visible symptom.
fn rad_quant_never_loses_fragments() {
    let tmp = tempfile::tempdir().unwrap();
    let rad = tmp.path().join("maps.rad");
    let names = ["t0", "t1", "t2"];
    let lens = [1000u32, 1000, 1000];

    // 487 fragments across enough chunks that the reader has real work to hand
    // out (write_rad flushes one chunk per 64 records).
    let mut frags: Vec<Vec<(u32, Option<u16>, i32)>> = Vec::new();
    for (tid, n) in [(0u32, 137), (1, 211), (2, 89)] {
        for _ in 0..n {
            frags.push(vec![(tid, Some(200), 0)]);
        }
    }
    for _ in 0..50 {
        frags.push(vec![(0, Some(200), 0), (2, Some(200), 0)]);
    }
    let total = frags.len() as u64;
    write_rad(&rad, &names, &lens, RadProfile::Sketch, &frags, None);

    // A single worker is the worst case: with no sibling to drain the queue,
    // one lost race loses the whole file.
    for threads in [1usize, 1, 2, 8] {
        let out = tmp.path().join(format!("q{threads}_{}", rand_suffix()));
        std::fs::create_dir_all(&out).unwrap();
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        let res = pool.install(|| quantify_rad(&opts_for(&out), &rad).unwrap());

        assert_eq!(
            res.num_processed, total,
            "-p {threads}: reader dropped fragments ({} of {total})",
            res.num_processed
        );
        assert_eq!(res.num_mapped, total, "-p {threads}: mapped-count mismatch");
        assert!(
            res.num_eq_classes > 0,
            "-p {threads}: no equivalence classes were built"
        );
        let sum: f64 = res.counts.iter().sum();
        assert!(
            (sum - total as f64).abs() < 1e-6,
            "-p {threads}: mass not conserved: {sum} vs {total}"
        );
    }
}

/// Distinct output directories across loop iterations that reuse a thread count.
fn rand_suffix() -> u64 {
    use std::time::{SystemTime, UNIX_EPOCH};
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos() as u64
}

#[test]
/// The headline guarantee of the RAD path: the answer cannot depend on how many
/// threads happened to be used.
fn thread_count_invariant() {
    let tmp = tempfile::tempdir().unwrap();
    let rad = tmp.path().join("maps.rad");
    let names = ["t0", "t1", "t2"];
    let lens = [1000u32, 1000, 1000];
    let mut frags: Vec<Vec<(u32, Option<u16>, i32)>> = Vec::new();
    for (tid, n) in [(0u32, 137), (1, 211), (2, 89)] {
        for _ in 0..n {
            frags.push(vec![(tid, Some(200), 0)]);
        }
    }
    // some ambiguous ones for good measure
    for _ in 0..50 {
        frags.push(vec![(0, Some(200), 0), (2, Some(200), 0)]);
    }
    write_rad(&rad, &names, &lens, RadProfile::Sketch, &frags, None);

    let run = |threads: usize| -> Vec<f64> {
        let out = tmp.path().join(format!("quant_{threads}"));
        std::fs::create_dir_all(&out).unwrap();
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        pool.install(|| quantify_rad(&opts_for(&out), &rad).unwrap().counts)
    };

    let c1 = run(1);
    let c8 = run(8);
    let c1b = run(1);
    assert_eq!(c1.len(), c8.len());
    // Quantification is fully deterministic: the FLD is fixed before eq-class
    // assembly (no online phase, no RNG), collapsed classes are order-independent,
    // and EM runs from a uniform start — so the result is bit-identical regardless
    // of thread count, and reproducible run-to-run.
    assert_eq!(
        c1, c8,
        "result differs across thread counts (non-deterministic)"
    );
    assert_eq!(
        c1, c1b,
        "result differs across repeated single-threaded runs"
    );
}

/// A salmon RAD that bakes its FLD into the header must use that FLD, not derive
/// one. The records imply fragment length 200, but we bake a distribution peaked
/// at 500; if the baked FLD is used, the reported mean is ~500 (not ~200).
#[test]
/// A baked fragment-length distribution must be taken verbatim, or a requant
/// would quietly differ from the run that produced the file.
fn baked_fld_is_used_not_derived() {
    let tmp = tempfile::tempdir().unwrap();
    let rad = tmp.path().join("maps.rad");
    let names = ["t0", "t1"];
    let lens = [4000u32, 4000];

    let mut frags: Vec<Vec<(u32, Option<u16>, i32)>> = Vec::new();
    for (tid, n) in [(0u32, 200), (1, 300)] {
        for _ in 0..n {
            frags.push(vec![(tid, Some(200), 0)]); // implied length 200
        }
    }
    // Baked log-PMF: an (unnormalized) Gaussian centered at 500 (from_log_pmf
    // normalizes). Reserved slot length is 1024 (see `write_rad`).
    let baked: Vec<f64> = (0..1024)
        .map(|l| {
            let d = (l as f64 - 500.0) / 20.0;
            -0.5 * d * d
        })
        .collect();
    write_rad(
        &rad,
        &names,
        &lens,
        RadProfile::Sketch,
        &frags,
        Some(&baked),
    );

    let out = tmp.path().join("quant");
    std::fs::create_dir_all(&out).unwrap();
    let res = quantify_rad(&opts_for(&out), &rad).expect("quantify_rad");

    // counts are still exact (FLD doesn't affect unique mappings)…
    assert!(
        (res.counts[0] - 200.0).abs() < 1.0,
        "t0 = {}",
        res.counts[0]
    );
    assert!(
        (res.counts[1] - 300.0).abs() < 1.0,
        "t1 = {}",
        res.counts[1]
    );
    // …but the FLD must be the baked one (mean ~500), proving it was not derived
    // from the records' implied length of 200.
    assert!(
        (res.frag_len_mean - 500.0).abs() < 20.0,
        "expected baked FLD mean ~500, got {}",
        res.frag_len_mean
    );
}

/// Regression for #1062: with a baked FLD in play, `--fldMean`/`--fldSD` are not
/// consulted at all, so the run must at least say where its distribution came
/// from.
#[test]
/// And the metadata must say so, so a user can tell where the distribution came
/// from.
fn baked_fld_is_reported_as_the_frag_length_source() {
    let tmp = tempfile::tempdir().unwrap();
    let rad = tmp.path().join("maps.rad");
    let out = tmp.path().join("quant");
    std::fs::create_dir_all(&out).unwrap();
    write_baked_fld_rad(&rad);

    let res = quantify_rad(&opts_for(&out), &rad).expect("quantify_rad");
    assert_eq!(res.frag_len_source.as_str(), "rad_baked");
}

/// `--fldPolicy derive` must ignore the baked distribution and rebuild one from
/// the RAD's own fragment lengths. The fixture bakes a mean of ~500 while every
/// record implies 200, so the two are impossible to confuse.
#[test]
/// `--fldPolicy derive` must genuinely re-derive rather than falling back to the
/// baked values.
fn fld_policy_derive_ignores_the_baked_distribution() {
    let tmp = tempfile::tempdir().unwrap();
    let rad = tmp.path().join("maps.rad");
    let out = tmp.path().join("quant");
    std::fs::create_dir_all(&out).unwrap();
    write_baked_fld_rad(&rad);

    let mut opts = opts_for(&out);
    opts.fld_policy = FldPolicy::Derive;
    let res = quantify_rad(&opts, &rad).expect("quantify_rad");

    assert_eq!(res.frag_len_source.as_str(), "rad_derived");
    assert!(
        (res.frag_len_mean - 200.0).abs() < 20.0,
        "expected the derived FLD (~200), got {}",
        res.frag_len_mean
    );
    // Counts stay exact: the FLD does not affect uniquely-mapped fragments.
    assert!(
        (res.counts[0] - 200.0).abs() < 1.0,
        "t0 = {}",
        res.counts[0]
    );
    assert!(
        (res.counts[1] - 300.0).abs() < 1.0,
        "t1 = {}",
        res.counts[1]
    );
}

/// `--fldPolicy prior` must ignore both the baked distribution (mean ~500) and
/// the records' implied lengths (200), leaving `--fldMean`/`--fldSD` in sole
/// control — the setting that makes a fragment-length sensitivity analysis
/// actually perturb the model.
#[test]
/// `prior` must hand control to --fldMean/--fldSD alone.
fn fld_policy_prior_puts_fld_args_in_control() {
    let tmp = tempfile::tempdir().unwrap();
    let rad = tmp.path().join("maps.rad");
    let out = tmp.path().join("quant");
    std::fs::create_dir_all(&out).unwrap();
    write_baked_fld_rad(&rad);

    let mut opts = opts_for(&out);
    opts.fld_policy = FldPolicy::Prior;
    opts.fld_mean = 350.0;
    opts.fld_sd = 15.0;
    let res = quantify_rad(&opts, &rad).expect("quantify_rad");

    assert_eq!(res.frag_len_source.as_str(), "prior");
    assert!(
        (res.frag_len_mean - 350.0).abs() < 1.0,
        "expected the supplied prior (350), got {}",
        res.frag_len_mean
    );
}

/// The point of the issue: under the default policy, two runs that differ only
/// in their fragment-length prior are indistinguishable; under `prior` they are
/// not. Effective lengths are what the FLD actually drives, so compare those.
#[test]
/// And those arguments must actually matter under that policy — otherwise the
/// mode would appear to work while doing nothing.
fn fld_policy_prior_makes_different_priors_produce_different_results() {
    let run = |policy: FldPolicy, mean: f64| {
        let tmp = tempfile::tempdir().unwrap();
        let rad = tmp.path().join("maps.rad");
        let out = tmp.path().join("quant");
        std::fs::create_dir_all(&out).unwrap();
        write_baked_fld_rad(&rad);
        let mut opts = opts_for(&out);
        opts.fld_policy = policy;
        opts.fld_mean = mean;
        opts.fld_sd = 20.0;
        let res = quantify_rad(&opts, &rad).expect("quantify_rad");
        (res.eff_lengths.clone(), res.frag_len_mean)
    };

    let (baked_lo, _) = run(FldPolicy::Baked, 150.0);
    let (baked_hi, _) = run(FldPolicy::Baked, 450.0);
    assert_eq!(
        baked_lo, baked_hi,
        "baked policy must stay reproducible regardless of the prior"
    );

    let (prior_lo, mean_lo) = run(FldPolicy::Prior, 150.0);
    let (prior_hi, mean_hi) = run(FldPolicy::Prior, 450.0);
    assert!((mean_lo - 150.0).abs() < 1.0, "{mean_lo}");
    assert!((mean_hi - 450.0).abs() < 1.0, "{mean_hi}");
    assert_ne!(
        prior_lo, prior_hi,
        "prior policy must let --fldMean change the effective lengths"
    );
}

/// Under auto (`-l A`) library type, the derived FLD must be built from the
/// majority-orientation unique fragments. Here 300 opposite-strand pairs imply
/// length 200 and 60 same-strand pairs imply length 800; the FLD must reflect the
/// majority (≈200), excluding the minority bucket's 800s.
#[test]
/// Orientation is decided by majority: the minority bucket is mis-mapping noise
/// and must not pollute the fragment-length distribution.
fn auto_orientation_derives_fld_from_majority_bucket() {
    let tmp = tempfile::tempdir().unwrap();
    let rad = tmp.path().join("maps.rad");
    let names = ["t0", "t1"];
    let lens = [4000u32, 4000];

    let w = RadOutputWriter::create(
        &rad,
        &names,
        &lens,
        true,
        RadProfile::Sketch,
        1024,
        salmon_rad::ChunkCodec::None,
        &prov(),
    )
    .unwrap();
    let ctx: SalmonBulkContext = *w.context();
    let mut cb = FragmentChunkBuf::with_capacity(4096);
    let mut emit = |_i: usize, tid: u32, fl: u16, is_fw: bool, mate_fw: bool| {
        let rec = salmon_rad::SalmonBulkRecord::new(
            frag_map_type::MAPPED_PAIR,
            vec![RadHit {
                tid,
                is_fw,
                mate_fw,
                pos: 0,
                frag_len: fl,
                score: 0,
            }],
        );
        cb.write(&rec, &ctx).unwrap();
    };
    let mut i = 0;
    for _ in 0..300 {
        emit(i, 0, 200, true, false); // opposite-strand (inward), length 200
        i += 1;
    }
    for _ in 0..60 {
        emit(i, 1, 800, true, true); // same-strand, length 800 (minority)
        i += 1;
    }
    w.append_chunk_bytes(&cb.take_bytes().unwrap()).unwrap();
    w.finalize().unwrap();

    // `-l A`: the format is auto-detected from the RAD (order-independent tally of
    // unique pairs). All 360 read-1s are forward and 300/360 are inward, so this
    // resolves to ISF — and the 60 same-strand (non-inward) pairs are then dropped
    // as strand-incompatible (default `incompat_prior` = 0 ⇒ ignore), leaving 300
    // mapped. The FLD is still derived only from the majority (opposite-strand)
    // bucket, so its mean is ~200, not contaminated by the minority's 800s.
    let mut o = AlignQuantOptions::new(tmp.path().join("unused.rad"), tmp.path().to_path_buf());
    o.lib_type = "A".to_string();
    o.fld_mean = 250.0;
    o.fld_sd = 25.0;
    let out = tmp.path().join("quant");
    std::fs::create_dir_all(&out).unwrap();
    let res = quantify_rad(&o, &rad).expect("quantify_rad");

    assert_eq!(
        res.num_mapped, 300,
        "auto-detected ISF should filter the 60 same-strand minority"
    );
    assert!(
        (res.frag_len_mean - 200.0).abs() < 15.0,
        "expected FLD mean ~200 from the majority bucket, got {} (minority 800s leaked in?)",
        res.frag_len_mean
    );
}

/// A `library_format` baked into the header is authoritative under `-l A`: even
/// though these records (all inward, read-1 forward) would auto-detect as ISF, a
/// baked ISR makes every placement strand-incompatible, so all are filtered out.
/// Proves the reader consults the baked format rather than re-inferring.
#[test]
/// A baked library format was determined by the run that saw the reads, so under
/// `-l A` it must win over re-detection from the RAD.
fn baked_library_format_is_authoritative_under_auto() {
    let tmp = tempfile::tempdir().unwrap();
    let rad = tmp.path().join("maps.rad");
    let names = ["t0"];
    let lens = [4000u32];

    let mut w = RadOutputWriter::create(
        &rad,
        &names,
        &lens,
        true,
        RadProfile::Sketch,
        1024,
        salmon_rad::ChunkCodec::None,
        &prov(),
    )
    .unwrap();
    let ctx: SalmonBulkContext = *w.context();
    let mut cb = FragmentChunkBuf::with_capacity(4096);
    for _ in 0..200 {
        let rec = salmon_rad::SalmonBulkRecord::new(
            frag_map_type::MAPPED_PAIR,
            vec![RadHit {
                tid: 0,
                is_fw: true,
                mate_fw: false, // inward, read-1 forward ⇒ would auto-detect ISF
                pos: 0,
                frag_len: 200,
                score: 0,
            }],
        );
        cb.write(&rec, &ctx).unwrap();
    }
    w.append_chunk_bytes(&cb.take_bytes().unwrap()).unwrap();
    // Bake ISR (format_id 2) — read-1 antisense, the opposite strandedness of what
    // these records imply.
    w.set_library_format(2);
    w.finalize().unwrap();

    let mut o = AlignQuantOptions::new(tmp.path().join("unused.rad"), tmp.path().to_path_buf());
    o.lib_type = "A".to_string();
    let res = quantify_rad(&o, &rad).expect("quantify_rad");

    // Baked ISR ⇒ the ISF-observed inward/forward pairs are strand-incompatible and
    // (default `incompat_prior` = 0) dropped — so nothing maps.
    assert_eq!(
        res.num_mapped, 0,
        "baked ISR must be applied (rejecting ISF-observed pairs), not re-inferred as ISF"
    );
}

/// Positional bias with `--rad`: exercise the full bias path — collection pass,
/// abundance-aware posterior, expected-pos models, corrected effective lengths,
/// re-EM — and confirm it runs, conserves mass, and populates the positional bias
/// dump. (`-t` is required for any `--rad` bias model, as in reads mode.)
#[test]
/// Positional bias correction must work from RAD input, which carries positions
/// but no sequence.
fn pos_bias_rad_runs() {
    let tmp = tempfile::tempdir().unwrap();
    let rad = tmp.path().join("maps.rad");
    let names = ["t0", "t1", "t2"];
    let lens = [3000u32, 3000, 3000];
    // a transcriptome FASTA matching the references (sequence content is
    // irrelevant to positional bias, but `-t` is required and lengths must match).
    let fa = tmp.path().join("txome.fa");
    {
        use std::io::Write as _;
        let mut f = std::fs::File::create(&fa).unwrap();
        for (i, &l) in lens.iter().enumerate() {
            let seq: String = (0..l).map(|j| b"ACGT"[(j as usize) % 4] as char).collect();
            writeln!(f, ">{}\n{}", names[i], seq).unwrap();
        }
    }
    // unique proper pairs spread across each transcript (varied positions feed the
    // positional model); a few multimappers exercise the abundance-aware posterior.
    let mut frags: Vec<Vec<(u32, Option<u16>, i32)>> = Vec::new();
    for (tid, n) in [(0u32, 200), (1, 150), (2, 250)] {
        for _ in 0..n {
            frags.push(vec![(tid, Some(200), 0)]);
        }
    }
    for _ in 0..40 {
        frags.push(vec![(0, Some(200), 0), (2, Some(200), 0)]);
    }
    write_rad(&rad, &names, &lens, RadProfile::Sketch, &frags, None);

    let mut o = opts_for(&tmp.path().join("quant"));
    o.pos_bias = true;
    o.transcripts = Some(fa);
    std::fs::create_dir_all(&o.output_dir).unwrap();
    let res = quantify_rad(&o, &rad).expect("quantify_rad --posBias");

    // mass conserved, and the positional bias model was populated.
    let sum: f64 = res.counts.iter().sum();
    assert!((sum - 640.0).abs() < 1.0, "Σcounts = {sum}");
    assert_eq!(res.num_mapped, 640);
    assert!(
        !res.bias_dump.obs5_pos.is_empty() && !res.bias_dump.exp5_pos.is_empty(),
        "positional bias dump should be populated"
    );
}

fn bias_ref_seq_fixture(tmp: &Path) -> (std::path::PathBuf, Vec<Vec<u8>>) {
    let rad = tmp.join("maps.rad");
    let names = ["t0", "t1", "t2"];
    let lens = [3000u32, 3000, 3000];
    let mut frags: Vec<Vec<(u32, Option<u16>, i32)>> = Vec::new();
    for (tid, n) in [(0u32, 200), (1, 150), (2, 250)] {
        for _ in 0..n {
            frags.push(vec![(tid, Some(200), 0)]);
        }
    }
    write_rad(&rad, &names, &lens, RadProfile::Sketch, &frags, None);

    // Reference sequences supplied directly (transcript-id order), as the index
    // would provide them — no FASTA on disk.
    let ref_seqs: Vec<Vec<u8>> = lens
        .iter()
        .map(|&l| (0..l).map(|j| b"ACGT"[(j as usize) % 4]).collect())
        .collect();

    (rad, ref_seqs)
}

/// Bias correction can take its reference sequences directly via
/// `AlignQuantOptions::ref_seqs` instead of a `-t` FASTA — the path
/// `--deterministic` uses to hand phase-2 the index's own sequences. No
/// `transcripts` is set here.
#[test]
/// GC correction needs sequence, which the index can supply — so it must work
/// without the user re-passing a transcriptome FASTA.
fn gc_bias_rad_uses_ref_seqs_without_transcripts() {
    let tmp = tempfile::tempdir().unwrap();
    let (rad, ref_seqs) = bias_ref_seq_fixture(tmp.path());

    let mut o = opts_for(&tmp.path().join("quant"));
    o.gc_bias = true;
    o.ref_seqs = Some(salmon_core::RefSeqs::from_sequences(ref_seqs));
    assert!(o.transcripts.is_none(), "this path must not need -t");
    std::fs::create_dir_all(&o.output_dir).unwrap();
    let res = quantify_rad(&o, &rad).expect("quantify_rad --gcBias via ref_seqs");

    let sum: f64 = res.counts.iter().sum();
    assert!((sum - 600.0).abs() < 1.0, "Σcounts = {sum}");
    assert_eq!(res.num_mapped, 600);
    assert!(
        !res.bias_dump.obs_gc.is_empty() && !res.bias_dump.exp_gc.is_empty(),
        "GC bias dump should be populated from the supplied ref_seqs"
    );
}

/// The per-transcript correction must produce finite, bit-identical results
/// across thread counts.
#[test]
/// Thread invariance must survive with every bias model enabled, which is where
/// the order-independent accumulators are actually exercised.
fn all_bias_rad_is_thread_invariant() {
    let tmp = tempfile::tempdir().unwrap();
    let (rad, ref_seqs) = bias_ref_seq_fixture(tmp.path());

    let run = |label: &str, threads: usize| {
        let mut o = opts_for(&tmp.path().join(format!("quant_{label}")));
        o.seq_bias = true;
        o.gc_bias = true;
        o.pos_bias = true;
        o.ref_seqs = Some(salmon_core::RefSeqs::from_sequences(ref_seqs.clone()));
        assert!(o.transcripts.is_none(), "this path must not need -t");
        std::fs::create_dir_all(&o.output_dir).unwrap();
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| quantify_rad(&o, &rad).expect("quantify_rad via ref_seqs"))
    };

    let single = run("bias_1t", 1);
    let parallel = run("bias_8t", 8);
    assert!(
        single.eff_lengths.iter().all(|value| value.is_finite()),
        "single-thread effective lengths must be finite"
    );
    assert!(
        parallel.eff_lengths.iter().all(|value| value.is_finite()),
        "parallel effective lengths must be finite"
    );
    assert_eq!(single.eff_lengths, parallel.eff_lengths);
    assert_eq!(single.counts, parallel.counts);
    assert_eq!(single.tpm, parallel.tpm);
    assert_eq!(single.bias_dump.obs_gc, parallel.bias_dump.obs_gc);
    assert_eq!(single.bias_dump.exp_gc, parallel.bias_dump.exp_gc);
    assert_eq!(single.bias_dump.obs5_seq, parallel.bias_dump.obs5_seq);
    assert_eq!(single.bias_dump.obs3_seq, parallel.bias_dump.obs3_seq);
    assert_eq!(single.bias_dump.exp5_seq, parallel.bias_dump.exp5_seq);
    assert_eq!(single.bias_dump.exp3_seq, parallel.bias_dump.exp3_seq);
    assert_eq!(single.bias_dump.obs5_pos, parallel.bias_dump.obs5_pos);
    assert_eq!(single.bias_dump.obs3_pos, parallel.bias_dump.obs3_pos);
    assert_eq!(single.bias_dump.exp5_pos, parallel.bias_dump.exp5_pos);
    assert_eq!(single.bias_dump.exp3_pos, parallel.bias_dump.exp3_pos);

    let sum: f64 = parallel.counts.iter().sum();
    assert!((sum - 600.0).abs() < 1.0, "Σcounts = {sum}");
    assert_eq!(parallel.num_mapped, 600);
    assert!(
        !parallel.bias_dump.obs_gc.is_empty() && !parallel.bias_dump.exp_gc.is_empty(),
        "GC bias dump should be populated from the supplied ref_seqs"
    );
    assert!(
        !parallel.bias_dump.obs5_seq.is_empty()
            && !parallel.bias_dump.exp5_seq.is_empty()
            && !parallel.bias_dump.obs5_pos.is_empty()
            && !parallel.bias_dump.exp5_pos.is_empty(),
        "sequence and positional bias dumps should be populated"
    );
}

/// The `@PG` chain a BAM carries must parse back into the fields
/// `meta_info.json` reports, including a command line containing spaces and
/// colons — the characters most likely to break a naive split.
#[test]
fn source_program_lines_parse_back() {
    let lines = vec![
        "@PG\tID:bowtie2\tPN:bowtie2\tVN:2.4.5\tCL:bowtie2 -x idx -1 a.fq -2 b.fq --rg-id x:y"
            .to_string(),
        "@PG\tID:samtools\tPN:samtools\tPP:bowtie2\tVN:1.17\tDS:sorted".to_string(),
    ];
    let progs = salmon_align::parse_source_programs(&lines);
    assert_eq!(progs.len(), 2);

    assert_eq!(progs[0].id, "bowtie2");
    assert_eq!(progs[0].program_name.as_deref(), Some("bowtie2"));
    assert_eq!(progs[0].version.as_deref(), Some("2.4.5"));
    assert_eq!(
        progs[0].command_line.as_deref(),
        Some("bowtie2 -x idx -1 a.fq -2 b.fq --rg-id x:y"),
        "a command line's own colons must not truncate it"
    );
    assert_eq!(progs[0].previous_id, None);

    assert_eq!(progs[1].previous_id.as_deref(), Some("bowtie2"));
    assert_eq!(progs[1].description.as_deref(), Some("sorted"));

    // Non-@PG lines and unknown tags are ignored rather than rejected: provenance
    // that fails to parse is worse than provenance reporting what it understands.
    let mixed = vec![
        "@HD\tVN:1.6".to_string(),
        "@PG\tID:x\tZZ:unknown-tag".to_string(),
    ];
    let progs = salmon_align::parse_source_programs(&mixed);
    assert_eq!(progs.len(), 1);
    assert_eq!(progs[0].id, "x");
}
