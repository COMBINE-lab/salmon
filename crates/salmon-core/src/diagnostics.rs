//! Shared run diagnostics: bad-input red flags derived from end-of-run
//! aggregates (never per-fragment), surfaced to the log and to
//! `meta_info.json`'s `diagnostics` array by both the reads and alignment
//! quantification paths.

use serde::Serialize;

/// A structured run diagnostic. `code` is a stable, machine-readable identifier
/// (snake_case) so downstream tools can key on it (e.g. `low_mapping_rate`,
/// `no_fragments_mapped`, `library_type_mismatch`). `severity` is one of
/// `"warning"` | `"error"` | `"info"`.
#[derive(Debug, Clone, Serialize)]
pub struct Diagnostic {
    pub code: String,
    pub severity: String,
    pub message: String,
}

impl Diagnostic {
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
pub fn peak_rss_kb() -> u64 {
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| {
            s.lines().find_map(|l| {
                l.strip_prefix("VmHWM:")
                    .and_then(|v| v.split_whitespace().next())
                    .and_then(|n| n.parse::<u64>().ok())
            })
        })
        .unwrap_or(0)
}
