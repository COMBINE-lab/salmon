//! What a `--writeMappings` / `--writeBam` file has to contain.
//!
//! These tests are about the *output contract*, not about mapping: they build a
//! tiny index, map a handful of reads whose correct placement is obvious by
//! construction, and then assert on the records. What they defend is everything
//! a reader is entitled to assume and that no unit test can see, because it only
//! exists once the whole pipeline has run: that the header describes the file,
//! that every record's fields agree with each other, and that the tags are
//! actually there.
//!
//! The single most valuable assertion here is the one about CIGAR and read
//! length. A record whose CIGAR does not account for exactly as many read bases
//! as `SEQ` holds is rejected by every SAM reader in existence, and it is the
//! failure mode that any change to CIGAR derivation risks introducing.

use std::path::Path;

use salmon_index::{build, IndexBuildOptions};
use salmon_quant::mapping_record::ReadGroup;
use salmon_quant::{quantify, QuantOptions};

/// Options for one fixture run, beyond the mapping output itself.
#[derive(Default)]
struct Run {
    read_group: Option<&'static str>,
    /// Quantify the first mates on their own, exercising the single-end path.
    single_end: bool,
}

/// Build a two-transcript index, quantify the fixture reads against it with SAM
/// output on, and return the file's text.
fn run_with_mapping_output(dir: &Path, run: Run) -> String {
    let transcripts = dir.join("txp.fa");
    let (a, b) = fixture_transcripts();
    std::fs::write(&transcripts, format!(">txA desc here\n{a}\n>txB\n{b}\n")).unwrap();

    let (reads1, reads2) = fixture_reads(&a, &b);
    let r1 = dir.join("r1.fq");
    let r2 = dir.join("r2.fq");
    std::fs::write(&r1, reads1).unwrap();
    std::fs::write(&r2, reads2).unwrap();

    let index = dir.join("idx");
    let mut build_options = IndexBuildOptions::new(vec![transcripts], index.clone());
    build_options.threads = 1;
    build(&build_options).expect("building the index");

    let sam = dir.join("mappings.sam");
    let mut options = QuantOptions::new(index, dir.join("out"));
    if run.single_end {
        options.unmated = vec![r1];
    } else {
        options.mates1 = vec![r1];
        options.mates2 = vec![r2];
    }
    options.lib_type = "IU".to_string();
    options.num_threads = 1;
    options.write_mappings = Some(sam.clone());
    options.read_group = run.read_group.map(|spec| ReadGroup::parse(spec).unwrap());
    quantify(&options).expect("quantifying");

    std::fs::read_to_string(&sam).expect("reading the SAM output")
}

/// Two transcripts sharing no k-mers, so every read has one obvious home.
fn fixture_transcripts() -> (String, String) {
    // A deterministic, non-repeating sequence: a linear congruential walk over
    // the four bases, seeded differently per transcript so they do not share
    // 31-mers.
    fn walk(seed: u64, len: usize) -> String {
        let mut state = seed;
        (0..len)
            .map(|_| {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                b"ACGT"[((state >> 33) % 4) as usize] as char
            })
            .collect()
    }
    (walk(12345, 1200), walk(98765, 900))
}

/// Inward-facing pairs drawn from each transcript, one of them carrying a
/// substitution so the output has a record with a non-trivial `MD`, and one
/// homopolymer pair that maps nowhere so the unaligned path has something to
/// write.
fn fixture_reads(a: &str, b: &str) -> (String, String) {
    let complement = |base: u8| match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        _ => b'A',
    };
    let revcomp = |s: &str| -> String {
        s.bytes()
            .rev()
            .map(|base| complement(base) as char)
            .collect()
    };
    let mut r1 = String::new();
    let mut r2 = String::new();
    let mut emit = |name: &str, forward: &str, reverse: &str| {
        let quality = "I".repeat(forward.len());
        r1.push_str(&format!("@{name}/1\n{forward}\n+\n{quality}\n"));
        let quality = "I".repeat(reverse.len());
        r2.push_str(&format!("@{name}/2\n{reverse}\n+\n{quality}\n"));
    };
    for (index, start) in [40usize, 200, 480, 700].into_iter().enumerate() {
        emit(
            &format!("a{index}"),
            &a[start..start + 60],
            &revcomp(&a[start + 150..start + 210]),
        );
    }
    for (index, start) in [30usize, 300, 500].into_iter().enumerate() {
        emit(
            &format!("b{index}"),
            &b[start..start + 60],
            &revcomp(&b[start + 120..start + 180]),
        );
    }
    // One mismatched read, so MD and NM have something to say.
    let mut mutated: Vec<u8> = a[600..660].bytes().collect();
    mutated[25] = complement(mutated[25]);
    emit(
        "mismatch",
        std::str::from_utf8(&mutated).unwrap(),
        &revcomp(&a[750..810]),
    );
    // A fragment that cannot map anywhere in this fixture.
    emit("nowhere", &"A".repeat(60), &"C".repeat(60));
    // An orphan: read 1 comes from txA, read 2 from nowhere at all.
    emit("orphan", &a[820..880], &"AC".repeat(30));
    (r1, r2)
}

/// Split a SAM line into its mandatory fields plus the tags.
fn fields(line: &str) -> Vec<&str> {
    line.split('\t').collect()
}

fn tag<'a>(fields: &[&'a str], name: &str) -> Option<&'a str> {
    fields
        .iter()
        .skip(11)
        .find(|f| f.starts_with(name))
        .copied()
}

fn records(sam: &str) -> Vec<Vec<&str>> {
    sam.lines()
        .filter(|line| !line.starts_with('@'))
        .map(fields)
        .collect()
}

/// How many read bases a CIGAR accounts for. `M`, `I` and `S` consume the read;
/// `D` and `N` do not.
fn cigar_read_length(cigar: &str) -> usize {
    let mut total = 0;
    let mut digits = String::new();
    for c in cigar.chars() {
        if c.is_ascii_digit() {
            digits.push(c);
        } else {
            let len: usize = digits.parse().unwrap_or(0);
            digits.clear();
            if matches!(c, 'M' | 'I' | 'S' | '=' | 'X') {
                total += len;
            }
        }
    }
    total
}

/// The header has to describe the file: its version, its ordering, one `@SQ` per
/// reference with a checksum, and the program that wrote it.
#[test]
fn the_header_describes_the_file() {
    let dir = tempfile::tempdir().unwrap();
    let sam = run_with_mapping_output(dir.path(), Run::default());
    let mut lines = sam.lines();

    let hd = lines.next().expect("a header");
    assert_eq!(hd, "@HD\tVN:1.6\tSO:unsorted\tGO:query");

    let sq: Vec<&str> = sam.lines().filter(|l| l.starts_with("@SQ")).collect();
    assert_eq!(sq.len(), 2, "one @SQ per transcript");
    for line in &sq {
        assert!(line.contains("\tLN:"), "@SQ needs a length: {line}");
        assert!(line.contains("\tM5:"), "@SQ needs a checksum: {line}");
        assert!(line.contains("\tUR:file://"), "@SQ needs a URI: {line}");
    }
    // The reference name stops at the first space: the rest of a FASTA header is
    // a description, not part of the name.
    assert!(
        sq.iter().any(|l| l.contains("\tSN:txA\t")),
        "reference name should not include the FASTA description: {sq:?}"
    );

    assert!(
        sam.lines().any(|l| l.starts_with("@PG\tID:salmon")),
        "the header must record what wrote the file"
    );
    assert!(
        sam.lines().any(|l| l.starts_with("@CO")),
        "the ordering contract should be stated in the file"
    );
}

/// Every record has to be internally consistent. This is the check a SAM reader
/// makes, so it is the check that decides whether the output is usable at all.
#[test]
fn every_record_is_internally_consistent() {
    let dir = tempfile::tempdir().unwrap();
    let sam = run_with_mapping_output(dir.path(), Run::default());
    let records = records(&sam);
    assert!(!records.is_empty(), "the run produced no records");

    for record in &records {
        assert!(record.len() >= 11, "too few fields: {record:?}");
        let name = record[0];
        let flags: u16 = record[1].parse().expect("FLAG is a number");
        let cigar = record[5];
        let sequence = record[9];
        let quality = record[10];

        assert!(
            !name.ends_with("/1") && !name.ends_with("/2"),
            "the mate suffix has to go, or a reader cannot pair the mates: {name}"
        );
        assert_eq!(
            cigar_read_length(cigar),
            sequence.len(),
            "CIGAR and SEQ disagree for {name}: {cigar} vs {} bases",
            sequence.len()
        );
        assert_eq!(
            quality.len(),
            sequence.len(),
            "QUAL and SEQ disagree for {name}"
        );
        assert!(
            quality.bytes().all(|q| (33..=126).contains(&q)),
            "QUAL must be printable for {name}: {quality}"
        );
        // Every mate here is part of a pair, so the pairing bits have to be set
        // and the mate reference named.
        assert!(flags & 0x1 != 0, "reads were paired: {name} has {flags}");
        assert!(
            record[6] == "=" || record[6] == "*",
            "both mates are on one transcript: {name} has RNEXT {}",
            record[6]
        );
    }
}

/// The tags are the point of the exercise: without them the file says nothing
/// beyond where the read landed.
#[test]
fn records_carry_the_tags_readers_expect() {
    let dir = tempfile::tempdir().unwrap();
    let sam = run_with_mapping_output(dir.path(), Run::default());
    let records = records(&sam);

    for record in &records {
        // An unmapped record has no placement, so none of the placement tags.
        if record[1].parse::<u16>().unwrap() & 0x4 != 0 {
            continue;
        }
        for name in ["NH:i:", "HI:i:", "AS:i:", "ZW:f:", "NM:i:", "MD:Z:"] {
            assert!(
                tag(record, name).is_some(),
                "{} is missing {name}: {record:?}",
                record[0]
            );
        }
        // A paired record knows its mate's shape without anyone sorting the
        // file, when there is a mate with a shape to know: an orphan's mate has
        // no CIGAR and no mapping quality to report.
        if record[1].parse::<u16>().unwrap() & 0x8 == 0 {
            assert!(tag(record, "MC:Z:").is_some(), "{} lacks MC", record[0]);
            assert!(tag(record, "MQ:i:").is_some(), "{} lacks MQ", record[0]);
        }
    }

    // MAPQ has to distinguish unique placements from ambiguous ones. These
    // transcripts share no k-mers, so every fragment is placed once.
    for record in records.iter().filter(|r| {
        let flags: u16 = r[1].parse().unwrap();
        flags & 0x4 == 0
    }) {
        assert_eq!(
            tag(record, "NH:i:"),
            Some("NH:i:1"),
            "fixture transcripts should not be ambiguous: {record:?}"
        );
        assert_eq!(record[4], "255", "a unique placement gets MAPQ 255");
    }

    // The mismatched read is the one record whose MD is not a bare run length.
    let mismatched = records
        .iter()
        .find(|r| r[0] == "mismatch" && r[1].parse::<u16>().unwrap() & 0x40 != 0)
        .expect("the mismatched read should be in the output");
    assert_eq!(tag(mismatched, "NM:i:"), Some("NM:i:1"));
    let md = tag(mismatched, "MD:Z:").unwrap();
    assert!(
        md.bytes().any(|b| b.is_ascii_alphabetic()),
        "MD should name the reference base at the mismatch: {md}"
    );
}

/// Single-end runs go through a separate path in the processor and a separate
/// branch of the record builder, and get their own kind of wrong: a stray
/// pairing bit, or a mate tag describing a mate that does not exist.
#[test]
fn single_end_records_claim_no_mate() {
    let dir = tempfile::tempdir().unwrap();
    let sam = run_with_mapping_output(
        dir.path(),
        Run {
            single_end: true,
            ..Run::default()
        },
    );
    let records = records(&sam);
    assert!(
        !records.is_empty(),
        "the single-end run produced no records"
    );
    for record in &records {
        let flags: u16 = record[1].parse().unwrap();
        assert_eq!(
            flags & 0x1,
            0,
            "a single-end read is not paired: {} has {flags}",
            record[0]
        );
        // Every bit that describes a mate is meaningless without one.
        assert_eq!(
            flags & (0x2 | 0x8 | 0x20 | 0x40 | 0x80),
            0,
            "{} carries pairing flags: {flags}",
            record[0]
        );
        assert_eq!(record[6], "*", "no mate reference: {record:?}");
        assert_eq!(record[7], "0", "no mate position: {record:?}");
        assert_eq!(record[8], "0", "no template length: {record:?}");
        assert!(tag(record, "MC:Z:").is_none(), "{} has MC", record[0]);
        assert!(tag(record, "MQ:i:").is_none(), "{} has MQ", record[0]);
        if flags & 0x4 == 0 {
            assert!(tag(record, "NM:i:").is_some(), "{} lacks NM", record[0]);
        }
    }
}

/// `AS` has to describe the alignment in the same record, not a different one.
///
/// It used to be the mapper's own score while the CIGAR, `NM` and `MD` came from
/// re-aligning the placement, so a record could claim a good score beside an
/// alignment full of edits. Anything that reads both, and every filter on `AS`
/// does, was being told two different stories.
#[test]
fn the_score_describes_the_alignment_in_its_own_record() {
    let dir = tempfile::tempdir().unwrap();
    let sam = run_with_mapping_output(dir.path(), Run::default());
    for record in records(&sam) {
        let (Some(score), Some(edits)) = (tag(&record, "AS:i:"), tag(&record, "NM:i:")) else {
            // A record without a base-level alignment carries neither, which is
            // the honest answer rather than a mismatched pair.
            assert_eq!(
                tag(&record, "NM:i:").is_none(),
                tag(&record, "MD:Z:").is_none(),
                "NM and MD go together: {record:?}"
            );
            continue;
        };
        let score: i32 = score.trim_start_matches("AS:i:").parse().unwrap();
        let edits: i32 = edits.trim_start_matches("NM:i:").parse().unwrap();
        let cigar = record[5];
        // Ungapped alignments are the ones whose score follows directly from the
        // edit count: matched bases earn 2, each mismatch swings 6.
        if !cigar.contains('I') && !cigar.contains('D') && !cigar.contains('S') {
            let aligned = cigar_read_length(cigar) as i32;
            assert_eq!(
                score,
                2 * aligned - 6 * edits,
                "AS and NM disagree for {}: {cigar}",
                record[0]
            );
        }
    }
}

/// A read group has to appear both in the header and on every record, or the
/// two disagree and a merged file loses track of where its reads came from.
#[test]
fn the_read_group_reaches_the_header_and_every_record() {
    let dir = tempfile::tempdir().unwrap();
    let sam = run_with_mapping_output(
        dir.path(),
        Run {
            read_group: Some(r"ID:s1\tSM:sample1\tPL:ILLUMINA"),
            ..Run::default()
        },
    );
    assert!(
        sam.lines()
            .any(|l| l == "@RG\tID:s1\tSM:sample1\tPL:ILLUMINA"),
        "the @RG line should be in the header"
    );
    for record in records(&sam) {
        assert_eq!(
            tag(&record, "RG:Z:"),
            Some("RG:Z:s1"),
            "every record should point at the read group: {record:?}"
        );
    }
}
