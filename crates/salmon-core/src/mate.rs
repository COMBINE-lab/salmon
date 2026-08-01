//! Mate/fragment status, mirroring pufferfish's `util::MateStatus`.
//!
//! # Background: mates and fragments
//!
//! In paired-end sequencing one physical DNA *fragment* is read from both ends,
//! producing two reads called *mates*. Ideally both mates place consistently on
//! the same transcript ("proper pair"), which pins down the fragment's length
//! and is the strongest evidence available. Often only one mate places, leaving
//! an *orphan*: still usable, but weaker evidence, since the fragment length is
//! then unknown. Single-end libraries have one read per fragment by design.
//!
//! This enum records which of those situations a given mapping represents, and
//! it is consulted throughout scoring, library-type compatibility, and output.

use serde::{Deserialize, Serialize};

/// Which part(s) of a fragment a mapping accounts for.
///
/// `Serialize`/`Deserialize` are derived so this can be written into JSON
/// metadata; `Copy` because it is one byte and gets passed around constantly.
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
    ///
    /// Note this is deliberately *not* true for [`Self::SingleEnd`]: a
    /// single-end read is complete evidence for its library, whereas an orphan
    /// is a pair with half its evidence missing.
    pub fn is_orphan(&self) -> bool {
        matches!(self, MateStatus::PairedEndLeft | MateStatus::PairedEndRight)
    }

    /// True for a properly paired mapping (both mates placed consistently).
    pub fn is_paired(&self) -> bool {
        matches!(self, MateStatus::PairedEndPaired)
    }
}
