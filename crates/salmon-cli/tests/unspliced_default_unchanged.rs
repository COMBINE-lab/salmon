//! Without `--unspliced`, `salmon index` and `salmon quant` must behave
//! exactly as before the option existed.
//!
//! The byte-for-byte comparison against the previous release needs both
//! binaries and lives in `scripts/compare_default_outputs.sh` (index, eight
//! quant configurations, every output file). This test pins what can be
//! checked from one build:
//!
//! * `--unspliced none` is the default, bit for bit: both spellings build
//!   byte-identical indices;
//! * a default index has exactly the pre-feature file set and `info.json`
//!   keys (no `unspliced` block, no `t2g_3col.tsv`);
//! * quant on it, on both inference paths, writes exactly the pre-feature file
//!   set and `meta_info.json` keys (no `quant.genes.usa.tsv`, no `unspliced`
//!   block).
//!
//! The key and file lists were read from the outputs of salmon at
//! COMBINE-lab/salmon@cf7c8ef0 (develop, before this option).

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::process::Command;

const SALMON: &str = env!("CARGO_BIN_EXE_salmon");

const INDEX_FILES: &[&str] = &[
    "duplicate_clusters.tsv",
    "index.ctab",
    "index.ectab",
    "index.refinfo",
    "index.ssi",
    "index.ssi.mphf",
    "index.tct",
    "index.tdct",
    "info.json",
    "refseq.bin",
    "refseq_offsets.json",
];

const INFO_KEYS: &[&str] = &[
    "k",
    "m",
    "canonical",
    "has_ec_table",
    "num_refs",
    "first_decoy_index",
    "num_decoys",
    "index_version",
    "sshash_format_version",
    "keep_duplicates",
    "salmon_version",
    "seq_hash",
    "name_hash",
    "seq_hash512",
    "name_hash512",
    "decoy_seq_hash",
    "decoy_name_hash",
];

const QUANT_FILES: &[&str] = &[
    "aux_info/ambig_info.tsv",
    "aux_info/expected_bias.gz",
    "aux_info/fld.gz",
    "aux_info/meta_info.json",
    "aux_info/observed_bias.gz",
    "aux_info/observed_bias_3p.gz",
    "cmd_info.json",
    "libParams/flenDist.txt",
    "lib_format_counts.json",
    "logs/salmon_quant.log",
    "quant.sf",
];

const META_KEYS: &[&str] = &[
    "salmon_version",
    "samp_type",
    "opt_type",
    "quant_errors",
    "num_libraries",
    "library_types",
    "frag_dist_length",
    "frag_length_mean",
    "frag_length_sd",
    "frag_length_source",
    "seq_bias_correct",
    "gc_bias_correct",
    "pos_bias_correct",
    "num_bias_bins",
    "mapping_type",
    "keep_duplicates",
    "index_seq_hash",
    "index_name_hash",
    "index_seq_hash512",
    "index_name_hash512",
    "index_decoy_seq_hash",
    "index_decoy_name_hash",
    "num_valid_targets",
    "num_decoy_targets",
    "num_eq_classes",
    "serialized_eq_classes",
    "eq_class_properties",
    "length_classes",
    "num_processed",
    "num_mapped",
    "num_dovetail_fragments",
    "num_fragments_filtered_vm",
    "num_alignments_below_threshold_for_mapped_fragments_vm",
    "percent_mapped",
    "num_decoy_fragments",
    "inference_truncated_mass",
    "num_bootstraps",
    "num_orphan",
    "range_factorization_bins",
    "meta_info_complete",
    "num_em_iterations",
    "em_converged",
    "detected_library_type",
    "inference_path",
    "total_time_seconds",
    "peak_rss_kb",
    "diagnostics",
    "call",
    "start_time",
    "end_time",
];

fn sequence(mut seed: u64, len: usize) -> String {
    (0..len)
        .map(|_| {
            seed = seed
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            "ACGT".as_bytes()[(seed >> 33) as usize & 3] as char
        })
        .collect()
}

/// Two transcripts plus a decoy, and inward paired reads from the transcripts.
fn fixture(dir: &Path) -> (PathBuf, PathBuf, PathBuf, PathBuf) {
    let a = sequence(3, 3000);
    let b = sequence(17, 3000);
    let decoy = sequence(41, 5000);
    let fasta = dir.join("gentrome.fa");
    std::fs::write(&fasta, format!(">txA\n{a}\n>txB\n{b}\n>chrD\n{decoy}\n")).unwrap();
    let decoys = dir.join("decoys.txt");
    std::fs::write(&decoys, "chrD\n").unwrap();
    let (mut r1, mut r2) = (String::new(), String::new());
    let qual = "I".repeat(100);
    for i in 0..400 {
        let source = if i % 3 == 0 { &a } else { &b };
        let start = (i * 13) % (source.len() - 300);
        let mate2: String = source[start + 200..start + 300]
            .chars()
            .rev()
            .map(|c| match c {
                'A' => 'T',
                'C' => 'G',
                'G' => 'C',
                _ => 'A',
            })
            .collect();
        r1.push_str(&format!(
            "@f{i}\n{}\n+\n{qual}\n",
            &source[start..start + 100]
        ));
        r2.push_str(&format!("@f{i}\n{mate2}\n+\n{qual}\n"));
    }
    let (p1, p2) = (dir.join("r1.fq"), dir.join("r2.fq"));
    std::fs::write(&p1, r1).unwrap();
    std::fs::write(&p2, r2).unwrap();
    (fasta, decoys, p1, p2)
}

fn salmon(args: &[&str]) {
    let o = Command::new(SALMON).args(args).output().unwrap();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
}

fn files(dir: &Path) -> BTreeSet<String> {
    fn walk(root: &Path, dir: &Path, out: &mut BTreeSet<String>) {
        for e in std::fs::read_dir(dir).unwrap() {
            let p = e.unwrap().path();
            if p.is_dir() {
                walk(root, &p, out);
            } else {
                let rel = p.strip_prefix(root).unwrap().to_string_lossy();
                out.insert(rel.replace('\\', "/"));
            }
        }
    }
    let mut out = BTreeSet::new();
    walk(dir, dir, &mut out);
    out
}

fn json_keys(path: &Path) -> Vec<String> {
    let v: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
    v.as_object().unwrap().keys().cloned().collect()
}

fn sorted(v: &[&str]) -> Vec<String> {
    let mut v: Vec<String> = v.iter().map(|s| s.to_string()).collect();
    v.sort();
    v
}

fn s(p: &Path) -> &str {
    p.to_str().unwrap()
}

#[test]
fn default_index_and_quant_are_untouched_by_the_unspliced_option() {
    let tmp = tempfile::tempdir().unwrap();
    let (fasta, decoys, r1, r2) = fixture(tmp.path());
    let (plain, none) = (tmp.path().join("plain"), tmp.path().join("none"));
    let base = ["index", "-t", s(&fasta), "-d", s(&decoys), "-p", "1", "-i"];
    salmon(&[&base[..], &[s(&plain)]].concat());
    salmon(&[&base[..], &[s(&none), "--unspliced", "none"]].concat());

    // `--unspliced none` is the default, bit for bit.
    let listed = files(&plain);
    assert_eq!(listed, files(&none));
    for f in &listed {
        assert!(
            std::fs::read(plain.join(f)).unwrap() == std::fs::read(none.join(f)).unwrap(),
            "{f} differs between the default and --unspliced none"
        );
    }
    // The pre-feature index layout and info.json keys.
    let expected: BTreeSet<String> = INDEX_FILES.iter().map(|f| f.to_string()).collect();
    assert_eq!(listed, expected);
    let mut keys = json_keys(&plain.join("info.json"));
    keys.sort();
    assert_eq!(keys, sorted(INFO_KEYS));

    // Quant on both inference paths: the pre-feature file set and keys.
    for (label, extra) in [("deterministic", None), ("online", Some("--online"))] {
        let out = tmp.path().join(label);
        let mut args = vec![
            "quant",
            "-i",
            s(&plain),
            "-l",
            "A",
            "-1",
            s(&r1),
            "-2",
            s(&r2),
            "-p",
            "2",
            "-o",
            s(&out),
        ];
        args.extend(extra);
        salmon(&args);
        let expected: BTreeSet<String> = QUANT_FILES.iter().map(|f| f.to_string()).collect();
        assert_eq!(files(&out), expected, "{label}");
        let mut keys = json_keys(&out.join("aux_info/meta_info.json"));
        keys.sort();
        assert_eq!(keys, sorted(META_KEYS), "{label}");
    }
}
