//! RAD-input quantification (`quantify_rad`): write a small salmon RAD with a
//! known mapping structure, quantify it, and check counts, mass conservation,
//! profile handling, and thread-count invariance.

use std::path::Path;

use salmon_align::{quantify_rad, AlignQuantOptions};
use salmon_rad::{
    frag_map_type, name_hash, FragmentChunkBuf, RadHit, RadOutputWriter, RadProfile,
    SalmonBulkContext, FRAG_LEN_UNPAIRED,
};

/// Write a salmon RAD file. Each fragment is a list of `(tid, frag_len)` hits;
/// `frag_len == None` marks an orphan/single-end placement.
fn write_rad(
    path: &Path,
    ref_names: &[&str],
    ref_lengths: &[u32],
    profile: RadProfile,
    fragments: &[Vec<(u32, Option<u16>, i32)>],
) {
    let w = RadOutputWriter::create(path, ref_names, ref_lengths, true, profile).unwrap();
    let ctx: SalmonBulkContext = *w.context();
    let mut cb = FragmentChunkBuf::with_capacity(4096);
    for (i, hits) in fragments.iter().enumerate() {
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
        let rec = salmon_rad::SalmonBulkRecord::new(
            frag_type,
            name_hash(format!("frag{i}").as_bytes()),
            rad_hits,
        );
        cb.write(&rec, &ctx).unwrap();
        // emit multiple chunks (real output flushes one chunk per minibatch)
        if cb.nrec() >= 64 {
            w.append_chunk_bytes(&cb.take_bytes()).unwrap();
        }
    }
    if cb.nrec() > 0 {
        w.append_chunk_bytes(&cb.take_bytes()).unwrap();
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

#[test]
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
    write_rad(&rad, &names, &lens, RadProfile::Sketch, &frags);

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
    write_rad(&rad, &names, &lens, RadProfile::Sketch, &frags);

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

/// Hand-write a piscem `map-bulk`-style RAD (read tag `frag_map_type`, alignment
/// triple `compressed_ori_ref`/`pos`/`frag_len`, no `frag_name_hash`) and confirm
/// salmon detects the piscem-bulk profile and quantifies it.
#[test]
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

#[test]
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
    write_rad(&rad, &names, &lens, RadProfile::Sketch, &frags);

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
    assert_eq!(c1.len(), c8.len());
    for (a, b) in c1.iter().zip(c8.iter()) {
        // collapsed equivalence classes are independent of chunk/worker order, so
        // the EM result must match regardless of thread count.
        assert!((a - b).abs() < 1e-6, "thread variance: {c1:?} vs {c8:?}");
    }
}
