//! Transcript-coverage profiles and degradation normalization.
//!
//! # The problem
//!
//! RNA in a sequencing library is not intact. It degrades, unevenly: longer
//! transcripts degrade faster than short ones, and how much any given sample
//! lost depends on how it was handled. A degraded transcript yields fewer
//! fragments than its abundance warrants, so its count is an underestimate —
//! and because the shortfall differs per transcript *and* per sample, no global
//! size factor removes it. A differential-expression test run on those counts
//! sees degradation differences as expression differences.
//!
//! DegNorm (Xiong *et al.*, *Genome Biology* 2019,
//! <https://nustatbioinfo.github.io/DegNorm/>) measures the effect directly
//! from the shape of each transcript's coverage curve and corrects the counts
//! by it. This crate implements that model on salmon's own data.
//!
//! # The two stages, and why there are two
//!
//! Degradation is only visible by comparison. One sample's coverage curve is
//! whatever it is; it takes *other* samples to establish what the curve should
//! have looked like and hence what is missing. So the work splits:
//!
//! * **During `salmon quant`** ([`coverage`]) — each fragment's placement is
//!   already in hand while the reads are streaming, so its contribution to the
//!   transcript's coverage curve is accumulated then and dumped to
//!   `aux_info/coverage.gz`. Doing it later would mean reading every read
//!   again.
//! * **Across the cohort** ([`cohort`], [`nmf`]) — `salmon degnorm` reads those
//!   dumps, fits the rank-one over-approximation per transcript, and writes the
//!   degradation indices and adjusted counts. The corrected counts are emitted
//!   both as matrices and as one `quant.sf` per sample, so an existing
//!   tximport/DESeq2 pipeline reads them without changing anything but the path
//!   it points at.
//!
//! # Differences from DegNorm
//!
//! Both tools fit the same model; they do not see the same data. DegNorm reads
//! genome BAMs and works per *gene*, on exon-union coverage with overlapping
//! exons removed, from primary alignments. salmon has no genome and no primary
//! alignment: it works per *transcript*, and a fragment that could have come
//! from several transcripts contributes to each in proportion to its posterior
//! probability. `docs/degnorm.md` spells the consequences out.

pub mod cohort;
pub mod coverage;
pub mod nmf;

pub use cohort::{run, write_adjusted_quants, write_tables, CohortOptions, CohortResult, Sample};
pub use coverage::{CoverageAccumulator, CoverageProfiles, DEFAULT_NUM_BINS};
pub use nmf::{fit, Fit, FitOptions};
