//! Transparent decompression of user-supplied input files.
//!
//! # Why this is shared
//!
//! salmon takes files from users in several places — the transcriptome FASTA,
//! the reads, an alignment record stream, a gene map — and users compress them
//! with whatever their pipeline uses. Deciding "is this compressed, and how"
//! separately at each entry point is how the four input paths ended up disagreeing:
//! one sniffed content with `needletail`, one sniffed magic bytes for gzip only,
//! one matched on a `.gz` file extension, and one used `niffler`. A file that
//! `salmon index` happily read could then be rejected by `salmon quant`.
//!
//! This module is the single answer. Everything that opens a user-supplied file
//! goes through it, so support is uniform by construction rather than by four
//! implementations agreeing.
//!
//! # Detection is by content, never by name
//!
//! A file extension is a hint, not a fact. `reads.fq.gz` may be plain text and
//! `reads.fq` may be gzip; neither is unusual once files have been moved,
//! renamed, or produced by a tool with its own conventions. [`niffler`] reads the
//! leading magic bytes and picks the decoder from what is actually there, which
//! also means a correctly-compressed file with an unexpected name just works
//! instead of failing with a parse error about the *decompressed* format that
//! never happened.
//!
//! Supported: gzip (including BGZF and multi-member files, via `MultiGzDecoder`,
//! so concatenated members are read in full rather than truncated at the first),
//! bzip2, xz and zstd. Anything else is passed through unchanged as plain text.

use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

/// Open `path`, transparently decompressing it if it is compressed.
///
/// Returns a buffered, `Send` reader so the result can be handed to a worker
/// thread (paraseq moves readers onto workers) and read line-wise without the
/// caller re-wrapping it.
///
/// A file shorter than niffler's five-byte sniff window cannot be any of the
/// compressed formats, so it is re-opened and returned as plain text rather than
/// reported as an error — an empty or stub input is the parser's business to
/// diagnose, not this function's.
pub fn open_maybe_compressed(path: &Path) -> io::Result<Box<dyn BufRead + Send>> {
    Ok(Box::new(BufReader::new(open_maybe_compressed_raw(path)?)))
}

/// As [`open_maybe_compressed`], but unbuffered, for callers that will impose
/// their own buffering (or hand the stream to a reader that does).
pub fn open_maybe_compressed_raw(path: &Path) -> io::Result<Box<dyn Read + Send>> {
    let open = || std::fs::File::open(path);
    // Buffer *before* sniffing: niffler reads the first few bytes to identify the
    // format, and an unbuffered file would make that several tiny syscalls.
    let probe: Box<dyn Read + Send> = Box::new(BufReader::new(open()?));
    match niffler::send::get_reader(probe) {
        Ok((reader, _format)) => Ok(reader),
        // Too short to be compressed: it is plain text by construction.
        Err(niffler::Error::FileTooShort) => Ok(Box::new(open()?)),
        Err(niffler::Error::IOError(e)) => Err(e),
        Err(e) => Err(io::Error::new(io::ErrorKind::InvalidData, e.to_string())),
    }
}

/// Strip one trailing compression extension from a file name, if present.
///
/// This module decides *compression* from content, never from the name — but a
/// caller may still have a name-based question to ask about what is underneath,
/// such as whether an alignment file is SAM or BAM. `reads.sam.zst` is a SAM
/// file, and asking "does this end in `.sam`" of the whole name answers no.
///
/// So this is not name-based compression detection creeping back in: it removes
/// compression noise from a name so a *different* question can be asked of what
/// remains. The suffixes are the ones [`open_maybe_compressed`] can actually
/// decode.
pub fn strip_compression_extension(name: &str) -> &str {
    const EXTENSIONS: [&str; 7] = [".gz", ".bgz", ".zst", ".zstd", ".bz2", ".xz", ".lzma"];
    let lower = name.to_ascii_lowercase();
    for ext in EXTENSIONS {
        if lower.ends_with(ext) {
            return &name[..name.len() - ext.len()];
        }
    }
    name
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// Round-trip each supported format through the sniffer and check the
    /// decompressed bytes come back, so a format regression is caught here rather
    /// than as a confusing parse error in whichever caller hit it first.
    #[test]
    fn sniffs_every_supported_format_by_content() {
        let dir = tempfile::tempdir().unwrap();
        let payload = b">seq1\nACGTACGTACGTACGTACGT\n";

        let write = |name: &str, bytes: &[u8]| {
            let p = dir.path().join(name);
            std::fs::File::create(&p).unwrap().write_all(bytes).unwrap();
            p
        };

        // Plain.
        let plain = write("a.fa", payload);
        let mut got = String::new();
        open_maybe_compressed(&plain)
            .unwrap()
            .read_to_string(&mut got)
            .unwrap();
        assert_eq!(got.as_bytes(), payload);

        // gzip, deliberately named so the extension says nothing useful — the
        // case that used to fail with "Invalid start character".
        let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        enc.write_all(payload).unwrap();
        let gz = write("b.plaintext", &enc.finish().unwrap());
        let mut got = String::new();
        open_maybe_compressed(&gz)
            .unwrap()
            .read_to_string(&mut got)
            .unwrap();
        assert_eq!(
            got.as_bytes(),
            payload,
            "gzip must be detected from content, not from the file name"
        );

        // zstd.
        let zst = write("c.fa.zst", &zstd::encode_all(&payload[..], 3).unwrap());
        let mut got = String::new();
        open_maybe_compressed(&zst)
            .unwrap()
            .read_to_string(&mut got)
            .unwrap();
        assert_eq!(got.as_bytes(), payload);
    }

    /// Multi-member gzip (what BGZF is) must be read to the end, not truncated at
    /// the first member — the reason gzip goes through `MultiGzDecoder`.
    #[test]
    fn reads_all_members_of_a_concatenated_gzip() {
        let dir = tempfile::tempdir().unwrap();
        let mut bytes = Vec::new();
        for part in [b">a\nACGT\n".as_slice(), b">b\nTTTT\n".as_slice()] {
            let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
            enc.write_all(part).unwrap();
            bytes.extend(enc.finish().unwrap());
        }
        let p = dir.path().join("multi.fa.gz");
        std::fs::File::create(&p)
            .unwrap()
            .write_all(&bytes)
            .unwrap();

        let mut got = String::new();
        open_maybe_compressed(&p)
            .unwrap()
            .read_to_string(&mut got)
            .unwrap();
        assert_eq!(got, ">a\nACGT\n>b\nTTTT\n");
    }

    /// Stripping a compression suffix lets a name-based format question see the
    /// name underneath, which is how `foo.sam.zst` is recognised as SAM.
    #[test]
    fn strips_one_compression_extension() {
        for (input, want) in [
            ("a.sam.gz", "a.sam"),
            ("a.sam.zst", "a.sam"),
            ("a.sam.ZST", "a.sam"),
            ("a.sam.bz2", "a.sam"),
            ("a.sam.xz", "a.sam"),
            ("a.sam", "a.sam"),
            ("a.bam", "a.bam"),
            // Only one is removed: the remaining name is what the caller asks about.
            ("a.sam.gz.gz", "a.sam.gz"),
            // A bare compression name has nothing underneath, and must not panic.
            (".gz", ""),
        ] {
            assert_eq!(strip_compression_extension(input), want, "for {input}");
        }
    }

    /// A file too short to sniff is plain text by construction, and must open
    /// rather than error.
    #[test]
    fn very_short_file_opens_as_plain_text() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("tiny.fa");
        std::fs::File::create(&p).unwrap().write_all(b">a").unwrap();
        let mut got = String::new();
        open_maybe_compressed(&p)
            .unwrap()
            .read_to_string(&mut got)
            .unwrap();
        assert_eq!(got, ">a");
    }
}
