//! Shared run diagnostics: bad-input red flags derived from end-of-run
//! aggregates (never per-fragment), surfaced to the log and to
//! `meta_info.json`'s `diagnostics` array by both the reads and alignment
//! quantification paths.
//!
//! # Why this exists
//!
//! salmon will happily produce a numerically valid `quant.sf` from reads that
//! have nothing to do with the reference: the numbers are just meaningless. The
//! commonest real-world failures are quiet ones — wrong organism, wrong
//! strandedness flag, untrimmed adapters — so the run ends by checking a few
//! end-of-run aggregates and saying plainly when the result looks untrustworthy.
//!
//! The checks deliberately use only totals, never per-fragment state, so they
//! cost nothing and cannot affect results.

use serde::Serialize;

/// A structured run diagnostic. `code` is a stable, machine-readable identifier
/// (snake_case) so downstream tools can key on it (e.g. `low_mapping_rate`,
/// `no_fragments_mapped`, `library_type_mismatch`). `severity` is one of
/// `"warning"` | `"error"` | `"info"`.
///
/// The code is what pipelines should match on; the `message` is prose and may be
/// reworded between releases.
#[derive(Debug, Clone, Serialize)]
pub struct Diagnostic {
    pub code: String,
    pub severity: String,
    pub message: String,
}

impl Diagnostic {
    /// Private constructor; the codes and severities are all chosen inside this
    /// module so they stay consistent.
    fn new(code: &str, severity: &str, message: String) -> Self {
        Self {
            code: code.to_string(),
            severity: severity.to_string(),
            message,
        }
    }
}

/// Build the standard bad-input diagnostics from already-computed aggregates.
/// Pure (no logging) so it is dependency-light; callers log the returned list.
///
/// - `num_processed` / `num_mapped`: fragment tallies.
/// - `auto_detect`: whether the library type was auto-detected (`-l A`); the
///   mismatch check only fires under an explicit `-l`.
/// - `requested_lib_type`: the resolved (used) library-type string.
/// - `detected_lib_type`: the format the detector observed (if any).
pub fn input_diagnostics(
    num_processed: u64,
    num_mapped: u64,
    auto_detect: bool,
    requested_lib_type: &str,
    detected_lib_type: Option<&str>,
) -> Vec<Diagnostic> {
    let mut d = Vec::new();
    // Ordered from most to least fundamental: no input at all, then input but no
    // mapping, then mapping but suspiciously little. Only one of the three can
    // apply, hence the if/else-if chain.
    if num_processed == 0 {
        d.push(Diagnostic::new(
            "no_input_fragments",
            "error",
            "no input fragments were processed — check that the input is non-empty and readable"
                .to_string(),
        ));
    } else if num_mapped == 0 {
        d.push(Diagnostic::new(
            "no_fragments_mapped",
            "error",
            "0 fragments mapped — the reference almost certainly does not match these reads (wrong reference/organism)"
                .to_string(),
        ));
    } else {
        let pct = 100.0 * num_mapped as f64 / num_processed as f64;
        // The thresholds are judgement calls, not statistics: a healthy bulk
        // RNA-seq run against a matching transcriptome typically maps 70-90%.
        // Below 30% something is usually wrong; below 10% it is almost certainly
        // the wrong reference.
        if pct < 10.0 {
            d.push(Diagnostic::new(
                "very_low_mapping_rate",
                "warning",
                format!("very low mapping rate ({pct:.1}%): the reference likely does not match these reads (wrong reference/organism or heavy contamination)"),
            ));
        } else if pct < 30.0 {
            d.push(Diagnostic::new(
                "low_mapping_rate",
                "warning",
                format!("low mapping rate ({pct:.1}%): verify the reference matches the sample and that adapter/quality trimming was applied"),
            ));
        }
    }
    // Strandedness sanity: an explicit `-l` disagreeing with the observed format.
    //
    // Under `-l A` there is nothing to disagree with — the detector's answer *is*
    // the setting — so the check is skipped, otherwise it would fire on itself.
    if !auto_detect {
        if let Some(det) = detected_lib_type {
            if det != requested_lib_type {
                d.push(Diagnostic::new(
                    "library_type_mismatch",
                    "warning",
                    format!("specified library type '{requested_lib_type}' disagrees with the observed format '{det}' — check the -l strandedness setting"),
                ));
            }
        }
    }
    d
}

/// Peak resident set size in KiB from `/proc/self/status` (`VmHWM`); 0 if
/// unavailable (non-Linux or parse failure). Read once at end of run.
///
/// "Peak RSS" is the high-water mark of physical memory the process ever held —
/// more useful than current usage for reporting, since the peak is what decides
/// whether a run fits on a given machine. Linux exposes it as the `VmHWM` line;
/// on other platforms the file simply does not exist and this returns 0.
pub fn peak_rss_kb() -> u64 {
    // Every step is fallible and every failure means the same thing ("no number
    // available"), so the chain uses `.ok()`/`and_then` to collapse them all
    // into a single `Option` and ends with a 0 default.
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| {
            s.lines().find_map(|l| {
                // Line looks like `VmHWM:\t  123456 kB`; take the first
                // whitespace-separated token after the prefix.
                l.strip_prefix("VmHWM:")
                    .and_then(|v| v.split_whitespace().next())
                    .and_then(|n| n.parse::<u64>().ok())
            })
        })
        .unwrap_or(0)
}

/// A `meta_info.json` field a complete description would carry, that this run
/// could not fill, and the upstream that could not supply it.
///
/// Recorded rather than silently omitted: a consumer reading a result long after
/// the run cannot otherwise tell "this run observed nothing" from "nobody could
/// tell this run". The `source` says which upstream to fix, so the record is
/// actionable and not merely an apology.
#[derive(Clone, Debug, serde::Serialize, PartialEq, Eq)]
pub struct MissingMetaField {
    /// the `meta_info.json` field that is absent or a placeholder
    pub field: &'static str,
    /// which upstream could not supply it
    pub source: MetaFieldSource,
    /// why it is unavailable, and where possible what would make it available
    pub reason: &'static str,
}

/// The upstream a [`MissingMetaField`] would have come from.
#[derive(Clone, Copy, Debug, serde::Serialize, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
pub enum MetaFieldSource {
    /// the salmon index the run was given
    Index,
    /// the RAD file the run quantified
    Rad,
    /// the BAM file the run quantified
    Bam,
}

impl MissingMetaField {
    /// Record a field the index could not supply.
    pub fn index(field: &'static str, reason: &'static str) -> Self {
        Self {
            field,
            source: MetaFieldSource::Index,
            reason,
        }
    }

    /// Record a field the input RAD could not supply.
    pub fn rad(field: &'static str, reason: &'static str) -> Self {
        Self {
            field,
            source: MetaFieldSource::Rad,
            reason,
        }
    }

    /// Record a field the input BAM could not supply.
    pub fn bam(field: &'static str, reason: &'static str) -> Self {
        Self {
            field,
            source: MetaFieldSource::Bam,
            reason,
        }
    }
}
