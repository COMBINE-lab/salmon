//! Access to reference sequences for alignment validation.
//!
//! The mapper needs the forward-strand bytes of a reference to align a read
//! against. piscem's index does not retain reference sequences, so they are
//! supplied through this trait (salmon-index persists and serves them). Keeping
//! it a trait decouples `salmon-map` from the concrete index/store.
//!
//! A "trait" in Rust is an interface: it lists methods a type must provide. Code
//! written against the trait works with any type that implements it, so the
//! mapper can be tested against a tiny in-memory stub and run against the real
//! on-disk index without changing a line.

/// Provides reference sequences (and decoy status) to the mapper.
pub trait RefProvider {
    /// Number of references.
    fn num_refs(&self) -> usize;

    /// Forward-strand sequence of reference `tid`.
    ///
    /// Returns a borrowed view rather than a copy: transcriptomes are hundreds
    /// of megabytes and this is called for every alignment.
    fn ref_view(&self, tid: u32) -> RefView<'_>;

    /// Whether reference `tid` is a decoy (genomic) sequence.
    ///
    /// Decoys are extra sequence (typically the genome) added to the index as a
    /// sink: a read that really came from an intron or an unannotated region
    /// will align better to the decoy than to any transcript, and is then
    /// discarded instead of being falsely credited to a transcript.
    ///
    /// This method has a body, so it is a *default*: an implementer that has no
    /// decoys can leave it out entirely.
    fn is_decoy(&self, _tid: u32) -> bool {
        false
    }
}

/// A stretch of reference bytes that two bits cannot hold (`N`, IUPAC codes,
/// lowercase): `len` copies of `byte` from position `start` of the concatenation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct OtherRun {
    pub start: u64,
    pub len: u32,
    pub byte: u8,
}

/// The base for each two-bit code; the code of a base is `(byte >> 1) & 3`.
const BASES: [u8; 4] = *b"ACTG";

/// The four bases held by each packed byte.
static UNPACKED: [[u8; 4]; 256] = {
    let mut t = [[0u8; 4]; 256];
    let mut i = 0;
    while i < 256 {
        t[i] = [
            BASES[i & 3],
            BASES[(i >> 2) & 3],
            BASES[(i >> 4) & 3],
            BASES[(i >> 6) & 3],
        ];
        i += 1;
    }
    t
};

#[inline]
fn is_acgt(b: u8) -> bool {
    (b == b'A') | (b == b'C') | (b == b'G') | (b == b'T')
}

/// Pack `bytes` into `out`, one byte per four bases, and append what two bits
/// cannot hold to `runs`. `start` is the position of `bytes[0]` in the
/// concatenation. `out` must hold `bytes.len().div_ceil(4)` bytes.
pub fn pack_bases(bytes: &[u8], out: &mut [u8], start: u64, runs: &mut Vec<OtherRun>) {
    let (chunks, rest) = bytes.as_chunks::<4>();
    for (c, o) in chunks.iter().zip(out.iter_mut()) {
        let w = u32::from_le_bytes(*c);
        let x = (w >> 1) & 0x0303_0303;
        let y = x | (x >> 6);
        *o = (y | (y >> 12)) as u8;
    }
    if !rest.is_empty() {
        let mut o = 0u8;
        for (i, &b) in rest.iter().enumerate() {
            o |= ((b >> 1) & 3) << (2 * i);
        }
        out[bytes.len() / 4] = o;
    }
    if bytes.iter().fold(true, |ok, &b| ok & is_acgt(b)) {
        return;
    }
    for (i, &b) in bytes.iter().enumerate() {
        if is_acgt(b) {
            continue;
        }
        let pos = start + i as u64;
        match runs.last_mut() {
            Some(r) if r.byte == b && r.start + u64::from(r.len) == pos && r.len < u32::MAX => {
                r.len += 1;
            }
            _ => runs.push(OtherRun {
                start: pos,
                len: 1,
                byte: b,
            }),
        }
    }
}

/// One reference's bases, stored as they are or at two bits per base.
#[derive(Clone, Copy)]
pub enum RefView<'a> {
    Bytes(&'a [u8]),
    Packed {
        /// the packed concatenation of all references
        bits: &'a [u8],
        /// where this reference starts in the concatenation
        start: u64,
        len: usize,
        /// the runs that overlap this reference
        other: &'a [OtherRun],
    },
}

impl<'a> RefView<'a> {
    /// A reference's view into a packed concatenation; `runs` covers all of it.
    #[inline]
    pub fn packed(bits: &'a [u8], start: u64, end: u64, runs: &'a [OtherRun]) -> Self {
        let other = if runs.is_empty() {
            runs
        } else {
            let lo = runs.partition_point(|r| r.start + u64::from(r.len) <= start);
            let hi = lo + runs[lo..].partition_point(|r| r.start < end);
            &runs[lo..hi]
        };
        Self::Packed {
            bits,
            start,
            len: (end - start) as usize,
            other,
        }
    }

    #[inline]
    pub fn len(&self) -> usize {
        match *self {
            Self::Bytes(b) => b.len(),
            Self::Packed { len, .. } => len,
        }
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The byte at position `i`, which must be in range.
    #[inline]
    pub fn base(&self, i: usize) -> u8 {
        match *self {
            Self::Bytes(b) => b[i],
            Self::Packed {
                bits,
                start,
                len,
                other,
            } => {
                debug_assert!(i < len);
                let abs = start + i as u64;
                let code = bits[(abs >> 2) as usize] >> ((abs & 3) * 2);
                if !other.is_empty() {
                    let k = other.partition_point(|r| r.start + u64::from(r.len) <= abs);
                    if let Some(r) = other.get(k) {
                        if r.start <= abs {
                            return r.byte;
                        }
                    }
                }
                BASES[(code & 3) as usize]
            }
        }
    }

    /// The bytes of `lo..hi`: borrowed where they are stored as they are, and
    /// otherwise decoded into `buf`.
    #[inline]
    pub fn slice<'b>(&self, lo: usize, hi: usize, buf: &'b mut Vec<u8>) -> &'b [u8]
    where
        'a: 'b,
    {
        match *self {
            Self::Bytes(b) => &b[lo..hi],
            Self::Packed {
                bits,
                start,
                len,
                other,
            } => {
                assert!(lo <= hi && hi <= len);
                decode(bits, start + lo as u64, start + hi as u64, other, buf)
            }
        }
    }

    /// How many leading bases of `query` match the reference from position `r`
    /// on, as [`Self::base`] compared without regard to case; `query` must fit.
    /// `None` when this view has no fast path and the bases must be compared one
    /// at a time.
    #[inline]
    pub fn match_right(&self, query: &[u8], r: usize) -> Option<usize> {
        let Self::Packed {
            bits,
            start,
            other: [],
            ..
        } = *self
        else {
            return None;
        };
        let mut abs = start as usize + r;
        let mut i = 0;
        while i < query.len() {
            let idx = abs >> 2;
            let mut w = match bits.get(idx..idx + 8) {
                Some(b) => u64::from_le_bytes(b.try_into().unwrap()),
                None => {
                    let mut t = [0u8; 8];
                    t[..bits.len() - idx].copy_from_slice(&bits[idx..]);
                    u64::from_le_bytes(t)
                }
            } >> ((abs & 3) * 2);
            let take = (32 - (abs & 3)).min(query.len() - i);
            for &q in &query[i..i + take] {
                if q & 0xDF != BASES[(w & 3) as usize] {
                    return Some(i);
                }
                w >>= 2;
                i += 1;
            }
            abs += take;
        }
        Some(i)
    }

    /// How many trailing bases of `query` match the reference just before
    /// position `r`; see [`Self::match_right`].
    #[inline]
    pub fn match_left(&self, query: &[u8], r: usize) -> Option<usize> {
        let Self::Packed {
            bits,
            start,
            other: [],
            ..
        } = *self
        else {
            return None;
        };
        // One past the base to compare next.
        let mut abs = start as usize + r;
        let mut i = 0;
        while i < query.len() {
            let idx = (abs - 1) >> 2;
            let lo = idx.saturating_sub(7);
            let mut t = [0u8; 8];
            t[7 - (idx - lo)..].copy_from_slice(&bits[lo..=idx]);
            // The next base on top.
            let p = (abs - 1) & 3;
            let mut w = u64::from_le_bytes(t) << ((3 - p) * 2);
            let take = ((idx - lo) * 4 + p + 1).min(query.len() - i);
            let end = query.len() - i;
            for &q in query[end - take..end].iter().rev() {
                if q & 0xDF != BASES[(w >> 62) as usize] {
                    return Some(i);
                }
                w <<= 2;
                i += 1;
            }
            abs -= take;
        }
        Some(i)
    }
}

/// The bytes of `lo..hi`, decoded into `out`, which they need not start.
fn decode<'b>(bits: &[u8], lo: u64, hi: u64, other: &[OtherRun], out: &'b mut Vec<u8>) -> &'b [u8] {
    // Whole packed bytes, trimmed to the bases asked for at the end.
    let skip = (lo & 3) as usize;
    let n = (hi - lo) as usize;
    let first = (lo >> 2) as usize;
    let last = (hi as usize).div_ceil(4);
    let want = (last - first) * 4;
    if out.len() < want {
        out.resize(want, 0);
    }
    let (chunks, _) = out[..want].as_chunks_mut::<4>();
    for (c, &b) in chunks.iter_mut().zip(&bits[first..last]) {
        *c = UNPACKED[b as usize];
    }
    let bases = &mut out[skip..skip + n];
    for r in other {
        let s = r.start.max(lo);
        let e = (r.start + u64::from(r.len)).min(hi);
        if s < e {
            bases[(s - lo) as usize..(e - lo) as usize].fill(r.byte);
        }
    }
    bases
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every window of a packed sequence reads back as the bytes it came from.
    #[test]
    fn packed_views_read_back_the_bytes() {
        let seq: Vec<u8> = b"ACGTTGCANNNNACGTacgtRYACGTAGGCTNACG".to_vec();
        for split in [4, 8, 12, 16] {
            let mut bits = vec![0u8; seq.len().div_ceil(4)];
            let mut runs = Vec::new();
            // Packed in two blocks, as the loader does.
            let (a, b) = seq.split_at(split);
            let (oa, ob) = bits.split_at_mut(split / 4);
            pack_bases(a, oa, 0, &mut runs);
            pack_bases(b, ob, split as u64, &mut runs);
            let mut buf = Vec::new();
            for start in 0..seq.len() {
                for end in start..=seq.len() {
                    let v = RefView::packed(&bits, start as u64, end as u64, &runs);
                    assert_eq!(v.len(), end - start);
                    for i in 0..v.len() {
                        assert_eq!(v.base(i), seq[start + i], "{start} {end} {i}");
                    }
                    for lo in 0..=v.len() {
                        for hi in lo..=v.len() {
                            assert_eq!(v.slice(lo, hi, &mut buf), &seq[start + lo..start + hi]);
                        }
                    }
                }
            }
        }
    }

    /// The word-at-a-time comparisons agree with comparing base by base.
    #[test]
    fn matches_agree_with_bases() {
        let mut x = 12345u64;
        let seq: Vec<u8> = (0..300)
            .map(|_| {
                x = x
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                b"ACGT"[(x >> 60) as usize & 3]
            })
            .collect();
        let mut bits = vec![0u8; seq.len().div_ceil(4)];
        pack_bases(&seq, &mut bits, 0, &mut Vec::new());
        let v = RefView::packed(&bits, 7, 290, &[]);
        for r in 0..v.len() {
            for take in [0, 1, 5, 33, 70, 200] {
                let n = take.min(v.len() - r);
                let mut q: Vec<u8> = seq[7 + r..7 + r + n].to_ascii_lowercase();
                for cut in [n, n / 2, n / 3] {
                    if cut < n {
                        q[cut] = b'N';
                    }
                    let want = q
                        .iter()
                        .enumerate()
                        .take_while(|&(i, &b)| b.to_ascii_uppercase() == v.base(r + i))
                        .count();
                    assert_eq!(v.match_right(&q, r), Some(want));
                }
                let n = take.min(r);
                let mut q: Vec<u8> = seq[7 + r - n..7 + r].to_vec();
                for cut in [n, n / 2, n / 3] {
                    if cut < n {
                        q[n - 1 - cut] = b'N';
                    }
                    let want = q
                        .iter()
                        .rev()
                        .enumerate()
                        .take_while(|&(i, &b)| b == v.base(r - 1 - i))
                        .count();
                    assert_eq!(v.match_left(&q, r), Some(want));
                }
            }
        }
    }
}
