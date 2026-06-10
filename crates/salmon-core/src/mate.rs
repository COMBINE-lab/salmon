//! Mate/fragment status, mirroring pufferfish's `util::MateStatus`.

use serde::{Deserialize, Serialize};

/// Which part(s) of a fragment a mapping accounts for.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum MateStatus {
    /// only the left/first mate of a pair mapped
    PairedEndLeft,
    /// only the right/second mate of a pair mapped
    PairedEndRight,
    /// both mates mapped as a proper pair
    PairedEndPaired,
    /// single-end read
    SingleEnd,
}

impl MateStatus {
    /// True when only one mate of a pair contributed (an orphan mapping).
    pub fn is_orphan(&self) -> bool {
        matches!(self, MateStatus::PairedEndLeft | MateStatus::PairedEndRight)
    }

    /// True for a properly paired mapping.
    pub fn is_paired(&self) -> bool {
        matches!(self, MateStatus::PairedEndPaired)
    }
}
