//! The decoder decision is made **per input file**, not once for the run.
//!
//! This matters for correctness, not just performance. A run may mix a regular
//! gzip file with a plain FASTQ, a zstd file, or a FIFO / process substitution
//! (`-1 <(zcat ...)`). The parallel decoder needs to seek, so it can only be
//! used on regular gzip files — but a single non-seekable input must not
//! downgrade the whole run, and must never be *probed*, because reading a byte
//! from a pipe consumes it and there is no second read.
//!
//! `open_input_pooled` decides per path: it checks `classify_input(path) ==
//! Regular` (metadata only — `std::fs::metadata`, so a FIFO is never touched)
//! and then whether the file is actually gzip. Anything failing either test
//! silently gets the serial decoder for that file alone. The run-level plan is
//! then reconciled against the decoder handles really obtained, rather than the
//! ones the decision hoped for.
//!
//! These tests pin that from salmon's side: the behaviour lives in piscem-rs,
//! and a regression there would otherwise show up here only as a slow run or,
//! worse, as reads silently consumed from a pipe.

use std::io::Write;
use std::path::PathBuf;

use piscem_rs::io::calibrate::DecoderPreference;
use salmon_quant::decode;
use thread_broker::EngagementPolicy;

/// Big enough that the budget clears any per-stream engagement threshold, so a
/// `false` below is always "could not", never "chose not to".
const BUDGET: usize = 64;

fn write_reads(dir: &std::path::Path, name: &str, gzip: bool) -> PathBuf {
    let mut fastq = Vec::new();
    for i in 0..2_000 {
        fastq.extend_from_slice(format!("@r{i}\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n").as_bytes());
    }
    let path = dir.join(name);
    if gzip {
        let mut enc =
            flate2::write::GzEncoder::new(std::fs::File::create(&path).unwrap(), Default::default());
        enc.write_all(&fastq).unwrap();
        enc.finish().unwrap();
    } else {
        std::fs::write(&path, &fastq).unwrap();
    }
    path
}

fn plan_for(paths: Vec<PathBuf>) -> decode::MappingPlan {
    // `Parallel` forces the decoder wherever it is possible at all, so what the
    // plan reports is purely "was this input decodable in parallel", with the
    // engagement policy taken out of the picture.
    decode::plan(
        &[paths],
        BUDGET,
        DecoderPreference::Parallel {
            workers_per_file: None,
        },
        EngagementPolicy::always(),
    )
    .expect("planning failed")
}

#[test]
fn both_gzip_inputs_get_the_parallel_decoder() {
    let dir = tempfile::tempdir().unwrap();
    let r1 = write_reads(dir.path(), "r1.fq.gz", true);
    let r2 = write_reads(dir.path(), "r2.fq.gz", true);
    let plan = plan_for(vec![r1, r2]);
    assert!(
        plan.parallel,
        "two regular gzip inputs should be decodable in parallel"
    );
    assert_eq!(plan.readers.len(), 2, "both inputs must still be opened");
}

#[test]
fn plain_inputs_fall_back_without_failing_the_run() {
    let dir = tempfile::tempdir().unwrap();
    let r1 = write_reads(dir.path(), "r1.fq", false);
    let r2 = write_reads(dir.path(), "r2.fq", false);
    let plan = plan_for(vec![r1, r2]);
    assert!(
        !plan.parallel,
        "uncompressed inputs yield no decoder handles, so the plan must \
         reconcile to serial rather than reserving slots for decoders that do \
         not exist"
    );
    assert_eq!(plan.readers.len(), 2, "both inputs must still be opened");
    assert_eq!(
        plan.decode_slots, 0,
        "slots reserved for a decoder that was never created are slots taken \
         away from mapping for nothing"
    );
}

/// The case the per-file rule exists for: one input can use the parallel
/// decoder and the other cannot.
#[test]
fn a_mixed_pair_downgrades_only_the_input_that_needs_it() {
    let dir = tempfile::tempdir().unwrap();
    let gz = write_reads(dir.path(), "r1.fq.gz", true);
    let plain = write_reads(dir.path(), "r2.fq", false);
    let plan = plan_for(vec![gz, plain]);

    assert_eq!(plan.readers.len(), 2, "both inputs must still be opened");
    assert!(
        plan.parallel,
        "one plain input must not cost the gzip input its parallel decoder: \
         the decision is per file, so the run keeps the decode concurrency it \
         can actually get"
    );
    // Sanity: the mapping side still holds the bulk of the budget.
    assert!(
        plan.map_threads > 0 && plan.map_threads + plan.decode_slots <= BUDGET,
        "split must stay within the budget: {} map + {} decode",
        plan.map_threads,
        plan.decode_slots
    );
}

/// A policy file overrides the per-mode default, and a typo in it is an error
/// rather than a silent no-op — a file that looks like it is configuring
/// something while doing nothing is worse than one that refuses to load.
#[test]
fn a_policy_file_overrides_the_mode_default() {
    let dir = tempfile::tempdir().unwrap();

    // No file: the mode decides, and the two modes must differ.
    let sa = decode::engagement_policy(false, None).unwrap();
    let sketch = decode::engagement_policy(true, None).unwrap();
    assert!(sa.min_threads_per_stream > sketch.min_threads_per_stream);

    // A file wins over the mode default, in both modes.
    let good = dir.path().join("policy.json");
    std::fs::write(&good, r#"{"parallel_decode": {"min_threads_per_stream": 3}}"#).unwrap();
    for sketch_mode in [true, false] {
        let p = decode::engagement_policy(sketch_mode, Some(&good)).unwrap();
        assert_eq!(
            p.min_threads_per_stream, 3,
            "an explicit policy must win over the per-mode default"
        );
    }

    let typo = dir.path().join("typo.json");
    std::fs::write(&typo, r#"{"parallel_decode": {"min_threads_per_input": 3}}"#).unwrap();
    assert!(
        decode::engagement_policy(false, Some(&typo)).is_err(),
        "a misspelled field must refuse to load rather than silently doing nothing"
    );
}
