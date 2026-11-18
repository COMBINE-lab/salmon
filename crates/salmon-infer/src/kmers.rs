//! GPU k-mer lookup for pseudoalignment (`--sketch`), after kallisto's GPU
//! port (Melsted et al., bioRxiv 2026.03.04.709526): a flat open-addressing
//! table from canonical k-mer to its unitig position, probed for every k-mer of
//! every read by one GPU thread per read.
//!
//! The kernel reproduces piscem's `SkippingStrategy::Strict` hit collection
//! (`HitSearcher::walk_safely_until`) exactly: it looks a k-mer up, emits a
//! raw hit, walks the read along that unitig while each next read k-mer is the
//! next unitig k-mer, and emits the walk's last matching k-mer as a second hit.
//! The CPU then runs piscem's target selection on those raw hits, so a read
//! maps to the same transcripts, at the same positions, as the CPU path.
//! [`KmerTable::hits`] is that kernel in Rust: the CPU reference the GPU is
//! tested against, and the definition of what it computes.
//!
//! Encoding (piscem's `ReadKmerIter`, sshash's SPSS): base `i` of a k-mer in
//! bits `2i..2i+2`, A=0 C=1 T=2 G=3 (`(b >> 1) & 3`, so lowercase works), the
//! complement `b ^ 2`; the canonical form is the smaller of a k-mer and its
//! reverse complement (`k` is odd, so they never tie).

use std::sync::atomic::{AtomicU32, AtomicU64, Ordering};

use rayon::prelude::*;
use std::borrow::Cow;

/// Raw hits a GPU thread can write per read; a read with more is redone by
/// the CPU (its count says so).
pub const HITS_CAP: usize = 32;

/// One raw hit: `[read_pos, unitig << 1 | fw_on_unitig, kmer_pos_in_unitig]`.
pub type Hit = [u32; 3];

const EMPTY: u64 = u64::MAX;

/// Canonical k-mer → `kmer_id << 1 | canonical_is_unitig_fw`, where
/// `kmer_id = kstart[unitig] + position`; linear probing at load ≤ 0.5.
pub struct KmerTable {
    k: u32,
    keys: Vec<u64>,
    vals: Vec<u32>,
    /// first k-mer id of each unitig, plus the total (`len = unitigs + 1`)
    kstart: Vec<u32>,
    /// the unitigs' bases, 2 bits each (base `b` in bits `2(b % 32)` of word
    /// `b / 32`), unitig `u` starting at base [`Self::base`]; a walk compares
    /// the read to it, as piscem does, instead of looking every k-mer up
    seq: Vec<u64>,
    /// `Some`: look k-mers up through piscem's own sshash dictionary
    /// structures instead of `keys`/`vals`/`kstart` (which are then empty)
    sshash: Option<Sshash>,
}

/// sshash's lookup structures, flattened for the GPU (option B): a k-mer's
/// minimizer selects a bucket, whose super-k-mer positions in the SPSS text
/// give the candidate k-mer starts, verified against `seq`. sshash's minimal
/// perfect hash minimizer → bucket is replaced by an exact open-addressing
/// table (an MPHF plus its presence check answers the same), and a heavy
/// bucket is scanned instead of going through the skew index (k-mers are
/// unique, so the scan finds the same one).
///
/// A table entry is `fp << 32 | bucket`, `fp` being the low half of the
/// minimizer's `fmix64` (whose high half picks the slot), so a probe reads one
/// word. The first entry with the query's `fp` along its probe sequence is the
/// query's own: the build refuses a table where another key's entry with the
/// same `fp` comes first. A minimizer absent from the dictionary may match an
/// `fp`, but no k-mer of the bucket it reaches can then match the read k-mer
/// (a dictionary k-mer's minimizer is a key), so the lookup still misses.
struct Sshash {
    m: u32,
    /// sshash's `DeterministicHasher` magic: `hash(v) = v * MIX ^ magic`
    magic: u64,
    /// each unitig's first base in `seq`, plus the total
    soff: Vec<u64>,
    /// bit `u` set: unitig `u` occurs in a reference (else its hits are dropped)
    keep: Vec<u32>,
    /// `fp << 32 | bucket` entries, [`EMPTY`] where free
    table: Vec<u64>,
    /// bucket `b`'s positions are `offs[bstart[b]..bstart[b + 1]]`
    bstart: Vec<u32>,
    /// super-k-mer minimizer positions in the SPSS text
    offs: Vec<u32>,
}

const MIX: u64 = 0x517c_c1b7_2722_0a95;

/// The `m`-mer at base `b` of `seq`.
#[inline]
fn word_at(seq: &[u64], b: u64, len: u32) -> u64 {
    let (w, sh) = ((b / 32) as usize, 2 * (b % 32));
    let mut x = seq[w] >> sh;
    if sh > 0 {
        x |= seq[w + 1] << (64 - sh);
    }
    x & ((1u64 << (2 * len)) - 1)
}

/// sshash's `compute_minimizer` for `k <= 31`: `(value, pos_in_kmer)`.
fn minimizer(fw: u64, rc: u64, k: u32, m: u32, magic: u64) -> (u64, u32) {
    let mask = (1u64 << (2 * m)) - 1;
    let kappa = |i: u32| ((fw >> (2 * i)) & mask).min((rc >> (2 * (k - m - i))) & mask);
    let hash = |v: u64| v.wrapping_mul(MIX) ^ magic;
    let mut min_value = kappa(0);
    let mut min_hash = hash(min_value);
    let (mut leftmost, mut rightmost, mut tie) = (0, 0, false);
    for i in 1..=k - m {
        let v = kappa(i);
        let h = hash(v);
        let (lt, eq) = (h < min_hash, h == min_hash);
        if lt {
            (min_hash, min_value, leftmost) = (h, v, i);
        }
        if lt || eq {
            rightmost = i;
        }
        tie = !lt && (tie || eq);
    }
    if !tie {
        return (min_value, leftmost);
    }
    // sshash's `resolve_tie`: the tied locus nearest the k-mer's centre
    let two_c = u64::from(k - m);
    let (mut best, mut lo, mut hi) = (u64::MAX, leftmost, leftmost);
    for i in leftmost..=rightmost {
        if hash(kappa(i)) == min_hash {
            let dist = (2 * u64::from(i)).abs_diff(two_c);
            if dist < best {
                (best, lo, hi) = (dist, i, i);
            } else if dist == best {
                hi = i;
            }
        }
    }
    (min_value, if lo == hi || rc >= fw { lo } else { hi })
}

#[inline]
fn fmix64(mut h: u64) -> u64 {
    h ^= h >> 33;
    h = h.wrapping_mul(0xff51_afd7_ed55_8ccd);
    h ^= h >> 33;
    h = h.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
    h ^ (h >> 33)
}

#[inline]
fn slot(key: u64, cap: u32) -> u32 {
    slot_fp(key, cap).0
}

/// The home slot of `key` and its fingerprint (see [`Sshash`]).
#[inline]
fn slot_fp(key: u64, cap: u32) -> (u32, u32) {
    let h = fmix64(key);
    ((((h >> 32) * u64::from(cap)) >> 32) as u32, h as u32)
}

/// Reverse complement of a `k`-mer word (piscem's `reverse_complement_word`).
#[inline]
pub fn revcomp(fw: u64, k: u32) -> u64 {
    let r = fw.reverse_bits();
    let r = ((r & 0x5555_5555_5555_5555) << 1) | ((r >> 1) & 0x5555_5555_5555_5555);
    (r >> (64 - 2 * k)) ^ (0xAAAA_AAAA_AAAA_AAAA >> (64 - 2 * k))
}

impl KmerTable {
    /// Index every k-mer of the unitigs `keep` accepts. `kstart` holds each
    /// unitig's first k-mer id and the total; `kmer(u, j)` is the `j`-th
    /// k-mer of unitig `u`, read in the unitig's own direction.
    ///
    /// `None` when `k` is even or above 31, the ids do not fit 31 bits, or a
    /// canonical k-mer occurs twice (never in a compacted de Bruijn graph).
    pub fn build(
        k: usize,
        kstart: Vec<u32>,
        keep: impl Fn(usize) -> bool + Sync,
        kmer: impl Fn(usize, u32) -> u64 + Sync,
    ) -> Option<Self> {
        let k = u32::try_from(k).ok().filter(|&k| k % 2 == 1 && k <= 31)?;
        let total = *kstart.last()?;
        if total >= 1 << 31 {
            return None;
        }
        let nu = kstart.len() - 1;
        let kept: u64 = (0..nu)
            .into_par_iter()
            .filter(|&u| keep(u))
            .map(|u| u64::from(kstart[u + 1] - kstart[u]))
            .sum();
        let cap = u32::try_from((2 * kept).max(1)).ok()?;
        let bases = u64::from(total) + nu as u64 * u64::from(k - 1);
        // one spare word, so a k-mer read never runs off the end
        let seq: Vec<AtomicU64> = (0..bases / 32 + 2).map(|_| AtomicU64::new(0)).collect();
        let keys: Vec<AtomicU64> = (0..cap).map(|_| AtomicU64::new(EMPTY)).collect();
        let vals: Vec<AtomicU32> = (0..cap).map(|_| AtomicU32::new(0)).collect();
        let unique = (0..nu).into_par_iter().filter(|&u| keep(u)).all(|u| {
            let base0 = u64::from(kstart[u]) + u as u64 * u64::from(k - 1);
            let put = |b: u64, x: u64| {
                seq[(b / 32) as usize].fetch_or(x << (2 * (b % 32)), Ordering::Relaxed);
            };
            (0..kstart[u + 1] - kstart[u]).all(|j| {
                let fw = kmer(u, j);
                if j == 0 {
                    (0..u64::from(k)).for_each(|i| put(base0 + i, (fw >> (2 * i)) & 3));
                } else {
                    put(base0 + u64::from(j + k - 1), fw >> (2 * (k - 1)));
                }
                let rc = revcomp(fw, k);
                let key = fw.min(rc);
                let mut s = slot(key, cap);
                loop {
                    match keys[s as usize].compare_exchange(
                        EMPTY,
                        key,
                        Ordering::Relaxed,
                        Ordering::Relaxed,
                    ) {
                        Ok(_) => break,
                        Err(other) if other == key => return false,
                        Err(_) => s = if s + 1 == cap { 0 } else { s + 1 },
                    }
                }
                let id = kstart[u] + j;
                vals[s as usize].store(id << 1 | u32::from(key == fw), Ordering::Relaxed);
                true
            })
        });
        // rayon's join is the happens-before edge for the relaxed stores above
        unique.then(|| Self {
            k,
            keys: keys.into_iter().map(AtomicU64::into_inner).collect(),
            vals: vals.into_iter().map(AtomicU32::into_inner).collect(),
            kstart,
            seq: seq.into_iter().map(AtomicU64::into_inner).collect(),
            sshash: None,
        })
    }

    /// Option B: look k-mers up through piscem's sshash dictionary. `spss` is
    /// its 2-bit text (base `i` in bits `2(i % 4)` of byte `i / 4`), `soff`
    /// each string's first base and the total, `keep(u)` whether unitig `u`
    /// occurs in a reference, `magic` the dictionary hasher's, and bucket `b`
    /// holds the super-k-mer positions `offs[bstart[b]..bstart[b + 1]]`.
    ///
    /// `None` when `k` is even or above 31, the text does not fit `u32`
    /// positions, or two buckets share a minimizer (never in sshash).
    #[allow(clippy::too_many_arguments)]
    pub fn build_sshash(
        k: usize,
        m: usize,
        magic: u64,
        spss: &[u8],
        soff: Vec<u64>,
        keep: impl Fn(usize) -> bool + Sync,
        bstart: Vec<u32>,
        offs: Vec<u32>,
    ) -> Option<Self> {
        let k = u32::try_from(k).ok().filter(|&k| k % 2 == 1 && k <= 31)?;
        let m = u32::try_from(m).ok().filter(|&m| m >= 1 && m <= k)?;
        if *soff.last()? >= 1 << 32 {
            return None;
        }
        // one spare word, so a k-mer read never runs off the end
        let mut seq = vec![0u64; spss.len().div_ceil(8) + 2];
        seq.par_iter_mut()
            .zip(spss.par_chunks(8))
            .for_each(|(w, chunk)| {
                let mut b = [0u8; 8];
                b[..chunk.len()].copy_from_slice(chunk);
                *w = u64::from_le_bytes(b);
            });
        let nu = soff.len() - 1;
        let mut bits = vec![0u32; nu.div_ceil(32).max(1)];
        let kept: Vec<bool> = (0..nu).into_par_iter().map(&keep).collect();
        for (u, _) in kept.iter().enumerate().filter(|(_, &x)| x) {
            bits[u / 32] |= 1 << (u % 32);
        }
        let nb = bstart.len() - 1;
        let cap = u32::try_from((2 * nb).max(1)).ok()?;
        let table: Vec<AtomicU64> = (0..cap)
            .into_par_iter()
            .map(|_| AtomicU64::new(EMPTY))
            .collect();
        // `b < nb < cap`, so no entry equals EMPTY
        let unique = (0..nb).into_par_iter().all(|b| {
            if bstart[b] == bstart[b + 1] {
                return true;
            }
            let mmer = word_at(&seq, u64::from(offs[bstart[b] as usize]), m);
            let (mut s, fp) = slot_fp(mmer.min(revcomp(mmer, m)), cap);
            let entry = u64::from(fp) << 32 | b as u64;
            loop {
                match table[s as usize].compare_exchange(
                    EMPTY,
                    entry,
                    Ordering::Relaxed,
                    Ordering::Relaxed,
                ) {
                    Ok(_) => return true,
                    // a duplicate key, or another key's `fp` ahead of this
                    // one's (insertion fills the first free slot, so every
                    // entry ahead on this probe sequence is already here)
                    Err(other) if (other >> 32) as u32 == fp => return false,
                    Err(_) => s = if s + 1 == cap { 0 } else { s + 1 },
                }
            }
        });
        unique.then(|| Self {
            k,
            keys: Vec::new(),
            vals: Vec::new(),
            kstart: Vec::new(),
            seq,
            sshash: Some(Sshash {
                m,
                magic,
                soff,
                keep: bits,
                table: table.into_iter().map(AtomicU64::into_inner).collect(),
                bstart,
                offs,
            }),
        })
    }

    /// Unitig `u`'s first base in `seq`.
    #[inline]
    fn base(&self, u: usize) -> u64 {
        match &self.sshash {
            Some(s) => s.soff[u],
            None => u64::from(self.kstart[u]) + u as u64 * u64::from(self.k - 1),
        }
    }

    /// `(unitig, k-mers in it, position, forward on it)` of the read k-mer.
    #[inline]
    fn locate(&self, fw: u64, rc: u64) -> Option<(usize, i64, i64, bool)> {
        let Some(s) = &self.sshash else {
            let (id, fw_on) = self.lookup(fw, rc)?;
            let u = self.unitig(id);
            let nk = i64::from(self.kstart[u + 1] - self.kstart[u]);
            return Some((u, nk, i64::from(id - self.kstart[u]), fw_on));
        };
        let k = self.k;
        let (v, p) = minimizer(fw, rc, k, s.m, s.magic);
        let cap = s.table.len() as u32;
        let (mut slot_, fp) = slot_fp(v, cap);
        let b = loop {
            match s.table[slot_ as usize] {
                EMPTY => return None,
                x if (x >> 32) as u32 == fp => break x as u32 as usize,
                _ => slot_ = if slot_ + 1 == cap { 0 } else { slot_ + 1 },
            }
        };
        for &pos in &s.offs[s.bstart[b] as usize..s.bstart[b + 1] as usize] {
            // sshash's `lookup_at_position`: the two loci that anchor the minimizer
            for back in [p, k - s.m - p] {
                let Some(kp) = u64::from(pos).checked_sub(u64::from(back)) else {
                    continue;
                };
                let stored = word_at(&self.seq, kp, k);
                let fw_on = if stored == fw {
                    true
                } else if stored == rc {
                    false
                } else {
                    continue;
                };
                let u = s.soff.partition_point(|&o| o <= kp) - 1;
                let nk = (s.soff[u + 1] - s.soff[u]) as i64 - i64::from(k) + 1;
                let cpos = (kp - s.soff[u]) as i64;
                if cpos >= nk {
                    continue;
                }
                // piscem skips hits on a unitig no reference contains
                if s.keep[u / 32] >> (u % 32) & 1 == 0 {
                    return None;
                }
                return Some((u, nk, cpos, fw_on));
            }
        }
        None
    }

    pub fn k(&self) -> u32 {
        self.k
    }

    /// `(kmer_id, k-mer is forward on its unitig)` of the read k-mer `fw`/`rc`.
    #[inline]
    fn lookup(&self, fw: u64, rc: u64) -> Option<(u32, bool)> {
        let key = fw.min(rc);
        let cap = self.keys.len() as u32;
        let mut s = slot(key, cap);
        loop {
            match self.keys[s as usize] {
                x if x == key => {
                    let v = self.vals[s as usize];
                    return Some((v >> 1, (fw == key) == (v & 1 == 1)));
                }
                EMPTY => return None,
                _ => s = if s + 1 == cap { 0 } else { s + 1 },
            }
        }
    }

    /// The unitig holding k-mer `id`.
    #[inline]
    fn unitig(&self, id: u32) -> usize {
        self.kstart.partition_point(|&s| s <= id) - 1
    }

    /// The kernel's buffers `(keys, vals, aux, words, prm)`. Option A: the
    /// table, `aux = kstart`, `words = seq`. Option B: the minimizer table,
    /// `aux = bstart ‖ offs ‖ keep`, `words = seq ‖ soff`. `prm`: k, table
    /// capacity, unitigs, hits per read, option B, m, magic (low, high), and
    /// where `offs`, `keep` and `soff` start.
    #[cfg_attr(not(gpu), allow(dead_code))]
    fn device(&self) -> (&[u64], &[u32], Cow<'_, [u32]>, Cow<'_, [u64]>, [u32; 11]) {
        let hcap = HITS_CAP as u32;
        let Some(s) = &self.sshash else {
            let nu = self.kstart.len() as u32 - 1;
            let prm = [
                self.k,
                self.keys.len() as u32,
                nu,
                hcap,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ];
            return (
                &self.keys,
                &self.vals,
                Cow::Borrowed(&self.kstart),
                Cow::Borrowed(&self.seq),
                prm,
            );
        };
        let aux = [&s.bstart[..], &s.offs, &s.keep].concat();
        let words = [&self.seq[..], &s.soff].concat();
        let prm = [
            self.k,
            s.table.len() as u32,
            s.soff.len() as u32 - 1,
            hcap,
            1,
            s.m,
            s.magic as u32,
            (s.magic >> 32) as u32,
            s.bstart.len() as u32,
            (s.bstart.len() + s.offs.len()) as u32,
            self.seq.len() as u32,
        ];
        // `vals` is unused here; one word, as CUDA refuses an empty upload
        (&s.table, &[0], Cow::Owned(aux), Cow::Owned(words), prm)
    }

    /// The raw hits of `read`, appended to `out`: the CPU reference of the GPU
    /// kernel ([`KMER_BODY`]), line for line.
    pub fn hits(&self, read: &[u8], out: &mut Vec<Hit>) {
        let k = self.k as i32;
        let mut it = Iter::new(read, k);
        let mut last_pos = -1;
        while !it.exhausted() {
            let p = it.pos as usize;
            if read[p..p + k as usize].iter().all(|&b| b == read[p]) {
                it.advance();
                continue;
            }
            let Some((u, nk, cpos, fw_on)) = self.locate(it.fw, it.rc) else {
                it.advance();
                continue;
            };
            let base = self.base(u);
            let start = it.pos;
            if start > last_pos {
                out.push([
                    start as u32,
                    (u as u32) << 1 | u32::from(fw_on),
                    cpos as u32,
                ]);
                last_pos = start;
            }
            let dir: i64 = if fw_on { 1 } else { -1 };
            let mut dist = if fw_on { nk - 1 - cpos } else { cpos };
            let mut ccur = cpos;
            let (mut lpos, mut lcpos, mut lfw, mut ended) = (start, cpos, fw_on, false);
            while !it.exhausted() && dist > 0 {
                let inc = i64::from(it.advance());
                dist -= inc;
                if dist < 0 {
                    break;
                }
                ccur += dir * inc;
                if it.exhausted() {
                    break;
                }
                let want = word_at(&self.seq, base + ccur as u64, self.k);
                match (it.fw == want, it.rc == want) {
                    (false, false) => break,
                    (fw2, _) => {
                        lfw = fw2;
                        lcpos += dir * inc;
                        lpos = it.pos;
                        ended = dist == 0;
                    }
                }
            }
            if lpos > last_pos {
                out.push([lpos as u32, (u as u32) << 1 | u32::from(lfw), lcpos as u32]);
                last_pos = lpos;
            }
            if ended || it.pos == start {
                it.advance();
            }
        }
    }
}

/// piscem's `ReadKmerIter`: the read's k-mer windows free of non-ACGT bytes.
struct Iter<'a> {
    read: &'a [u8],
    k: i32,
    pos: i32,
    last_inv: i32,
    fw: u64,
    rc: u64,
}

#[inline]
fn valid(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')
}

#[inline]
fn enc(b: u8) -> u64 {
    u64::from((b >> 1) & 3)
}

impl<'a> Iter<'a> {
    fn new(read: &'a [u8], k: i32) -> Self {
        let mut it = Self {
            read,
            k,
            pos: 0,
            last_inv: -1,
            fw: 0,
            rc: 0,
        };
        it.find_next_valid();
        if !it.exhausted() {
            it.build();
        }
        it
    }

    fn exhausted(&self) -> bool {
        self.pos > self.read.len() as i32 - self.k
    }

    fn build(&mut self) {
        let w = &self.read[self.pos as usize..(self.pos + self.k) as usize];
        self.fw = w
            .iter()
            .enumerate()
            .fold(0, |a, (i, &b)| a | enc(b) << (2 * i));
        self.rc = w.iter().fold(0, |a, &b| a << 2 | (enc(b) ^ 2));
    }

    fn find_next_valid(&mut self) {
        let len = self.read.len() as i32;
        loop {
            if self.pos + self.k > len {
                return;
            }
            if self.last_inv >= self.pos {
                self.pos = self.last_inv + 1;
                if self.pos + self.k > len {
                    return;
                }
            }
            let from = (self.last_inv + 1).max(self.pos);
            match (from..self.pos + self.k).find(|&i| !valid(self.read[i as usize])) {
                Some(i) => self.last_inv = i,
                None => return,
            }
        }
    }

    /// Move one window right (past any non-ACGT byte); the distance moved.
    fn advance(&mut self) -> i32 {
        let old = self.pos;
        self.pos += 1;
        let nr = self.pos + self.k - 1;
        if nr < self.read.len() as i32 && !valid(self.read[nr as usize]) {
            self.last_inv = nr;
        }
        if self.last_inv >= self.pos {
            self.pos = self.last_inv + 1;
            self.find_next_valid();
            if !self.exhausted() {
                self.build();
            }
        } else if !self.exhausted() {
            let b = enc(self.read[nr as usize]);
            self.fw = self.fw >> 2 | b << (2 * (self.k - 1));
            self.rc = (self.rc << 2 | (b ^ 2)) & ((1u64 << (2 * self.k)) - 1);
        }
        self.pos - old
    }
}

/// [`KmerTable::hits`] in C, shared by the Metal (MSL) and the CUDA/HIP
/// sources: each backend prepends `u8`/`u32`/`u64`/`i64`, `G` (the device
/// address space), `PRIV` (thread-private) and `FN` (function qualifiers),
/// and appends its entry point calling `map_mate`.
#[cfg_attr(not(gpu), allow(dead_code))]
pub(crate) const KMER_BODY: &str = r#"
FN u64 fmix64(u64 h) {
    h ^= h >> 33; h *= (u64)0xff51afd7ed558ccdUL;
    h ^= h >> 33; h *= (u64)0xc4ceb9fe1a85ec53UL;
    return h ^ (h >> 33);
}
FN bool valid_base(u8 b) {
    return b == 'A' || b == 'C' || b == 'G' || b == 'T' || b == 'a' || b == 'c' || b == 'g' || b == 't';
}
FN u64 enc(u8 b) { return (u64)((b >> 1) & 3); }

struct Iter { int pos; int last_inv; int len; int k; u64 fw; u64 rc; };

FN bool exhausted(PRIV const Iter& it) { return it.pos > it.len - it.k; }
FN void build_words(PRIV Iter& it, G const u8* r) {
    u64 fw = 0, rc = 0;
    for (int i = 0; i < it.k; ++i) {
        u64 b = enc(r[it.pos + i]);
        fw |= b << (2 * i);
        rc = (rc << 2) | (b ^ 2);
    }
    it.fw = fw; it.rc = rc;
}
FN void find_next_valid(PRIV Iter& it, G const u8* r) {
    for (;;) {
        if (it.pos + it.k > it.len) return;
        if (it.last_inv >= it.pos) {
            it.pos = it.last_inv + 1;
            if (it.pos + it.k > it.len) return;
        }
        int from = it.last_inv + 1 > it.pos ? it.last_inv + 1 : it.pos;
        bool bad = false;
        for (int i = from; i < it.pos + it.k; ++i) {
            if (!valid_base(r[i])) { it.last_inv = i; bad = true; break; }
        }
        if (!bad) return;
    }
}
FN int advance(PRIV Iter& it, G const u8* r) {
    int old = it.pos;
    it.pos += 1;
    int nr = it.pos + it.k - 1;
    if (nr < it.len && !valid_base(r[nr])) it.last_inv = nr;
    if (it.last_inv >= it.pos) {
        it.pos = it.last_inv + 1;
        find_next_valid(it, r);
        if (!exhausted(it)) build_words(it, r);
    } else if (!exhausted(it)) {
        u64 b = enc(r[nr]);
        it.fw = (it.fw >> 2) | (b << (2 * (it.k - 1)));
        it.rc = ((it.rc << 2) | (b ^ 2)) & ((((u64)1) << (2 * it.k)) - 1);
    }
    return it.pos - old;
}
// kmer id of the read k-mer, or -1; `fw_on`: forward on its unitig
FN i64 lookup(u64 fw, u64 rc, G const u64* keys, G const u32* vals, u32 cap, PRIV bool& fw_on) {
    u64 key = fw < rc ? fw : rc;
    u32 s = (u32)(((fmix64(key) >> 32) * (u64)cap) >> 32);
    for (;;) {
        u64 x = keys[s];
        if (x == key) {
            u32 v = vals[s];
            fw_on = (fw == key) == ((v & 1u) == 1u);
            return (i64)(v >> 1);
        }
        if (x == ~((u64)0)) return -1;
        s = s + 1 == cap ? 0 : s + 1;
    }
}
FN void emit(G u32* out, PRIV u32& nh, u32 hcap, int pos, u32 u, bool fw, i64 cpos) {
    if (nh < hcap) {
        out[3 * nh] = (u32)pos;
        out[3 * nh + 1] = (u << 1) | (fw ? 1u : 0u);
        out[3 * nh + 2] = (u32)cpos;
    }
    nh += 1;
}
// the `len`-mer at base `b` of the packed bases `w`
FN u64 word_at(G const u64* w, u64 b, int len) {
    u64 i = b / 32, sh = 2 * (b % 32);
    u64 x = w[i] >> sh;
    if (sh > 0) x |= w[i + 1] << (64 - sh);
    return x & ((((u64)1) << (2 * len)) - 1);
}
// option A: the k-mer table; the unitig, or -1
FN i64 locate_a(u64 fw, u64 rc, G const u64* keys, G const u32* vals, G const u32* kstart,
                G const u32* prm, PRIV i64& nk, PRIV i64& cpos, PRIV bool& fw_on) {
    i64 id = lookup(fw, rc, keys, vals, prm[1], fw_on);
    if (id < 0) return -1;
    u32 lo = 0, hi = prm[2];
    while (hi - lo > 1) {
        u32 mid = (lo + hi) / 2;
        if ((i64)kstart[mid] <= id) lo = mid; else hi = mid;
    }
    nk = (i64)(kstart[lo + 1] - kstart[lo]);
    cpos = id - (i64)kstart[lo];
    return (i64)lo;
}
FN u64 kappa(u64 fw, u64 rc, int k, int m, int i, u64 mask) {
    u64 a = (fw >> (2 * i)) & mask, b = (rc >> (2 * (k - m - i))) & mask;
    return a < b ? a : b;
}
FN u64 mini_hash(u64 v, u64 magic) { return (v * (u64)0x517cc1b727220a95UL) ^ magic; }
// sshash's compute_minimizer and resolve_tie: the minimizer, its position in `pos`
FN u64 minimizer(u64 fw, u64 rc, int k, int m, u64 magic, PRIV int& pos) {
    u64 mask = (((u64)1) << (2 * m)) - 1;
    u64 mv = kappa(fw, rc, k, m, 0, mask), mh = mini_hash(mv, magic);
    int left = 0, right = 0;
    bool tie = false;
    for (int i = 1; i <= k - m; ++i) {
        u64 v = kappa(fw, rc, k, m, i, mask), h = mini_hash(v, magic);
        bool lt = h < mh, eq = h == mh;
        if (lt) { mh = h; mv = v; left = i; }
        if (lt || eq) right = i;
        tie = !lt && (tie || eq);
    }
    pos = left;
    if (!tie) return mv;
    u64 best = ~((u64)0);
    int lo = left, hi = left;
    for (int i = left; i <= right; ++i) {
        if (mini_hash(kappa(fw, rc, k, m, i, mask), magic) == mh) {
            int d = 2 * i - (k - m);
            u64 dist = (u64)(d < 0 ? -d : d);
            if (dist < best) { best = dist; lo = i; hi = i; } else if (dist == best) hi = i;
        }
    }
    pos = (lo == hi || rc >= fw) ? lo : hi;
    return mv;
}
// option B: piscem's sshash dictionary (minimizer → bucket → super-k-mer
// positions, verified on the text); the unitig, or -1
FN i64 locate_b(u64 fw, u64 rc, G const u64* table, G const u32* aux,
                G const u64* words, G const u32* prm, PRIV i64& nk, PRIV i64& cpos, PRIV bool& fw_on) {
    int k = (int)prm[0], m = (int)prm[5];
    u32 cap = prm[1], nu = prm[2];
    G const u32* offs = aux + prm[8];
    G const u32* keep = aux + prm[9];
    G const u64* soff = words + prm[10];
    int p = 0;
    u64 v = minimizer(fw, rc, k, m, (u64)prm[6] | ((u64)prm[7] << 32), p);
    u64 h = fmix64(v);
    u32 s = (u32)(((h >> 32) * (u64)cap) >> 32), fp = (u32)h;
    u32 b;
    for (;;) {
        u64 x = table[s];
        if (x == ~((u64)0)) return -1;
        if ((u32)(x >> 32) == fp) { b = (u32)x; break; }
        s = s + 1 == cap ? 0 : s + 1;
    }
    for (u32 j = aux[b]; j < aux[b + 1]; ++j) {
        u64 pos = (u64)offs[j];
        for (int t = 0; t < 2; ++t) {
            u64 back = (u64)(t == 0 ? p : k - m - p);
            if (pos < back) continue;
            u64 kp = pos - back;
            u64 stored = word_at(words, kp, k);
            if (stored == fw) fw_on = true; else if (stored == rc) fw_on = false; else continue;
            u32 lo = 0, hi = nu;
            while (hi - lo > 1) {
                u32 mid = (lo + hi) / 2;
                if (soff[mid] <= kp) lo = mid; else hi = mid;
            }
            nk = (i64)(soff[lo + 1] - soff[lo]) - k + 1;
            cpos = (i64)(kp - soff[lo]);
            if (cpos >= nk) continue;
            if (((keep[lo / 32] >> (lo % 32)) & 1u) == 0) return -1;
            return (i64)lo;
        }
    }
    return -1;
}
// prm: see `KmerTable::device`
FN void map_mate(u32 i, G const u8* seq, G const u32* off, G const u64* keys, G const u32* vals,
                 G const u32* aux, G const u64* words, G const u32* prm, G u32* hits, G u32* counts) {
    G const u8* r = seq + off[i];
    int k = (int)prm[0];
    u32 hcap = prm[3];
    bool sshash = prm[4] != 0;
    G u32* out = hits + (u64)i * hcap * 3;
    u32 nh = 0;
    int last_pos = -1;
    Iter it;
    it.pos = 0; it.last_inv = -1; it.len = (int)(off[i + 1] - off[i]); it.k = k; it.fw = 0; it.rc = 0;
    find_next_valid(it, r);
    if (!exhausted(it)) build_words(it, r);
    while (!exhausted(it)) {
        bool homo = true;
        for (int j = 1; j < k; ++j) {
            if (r[it.pos + j] != r[it.pos]) { homo = false; break; }
        }
        if (homo) { advance(it, r); continue; }
        bool fw_on = false;
        i64 nk = 0, cpos = 0;
        i64 u = sshash ? locate_b(it.fw, it.rc, keys, aux, words, prm, nk, cpos, fw_on)
                       : locate_a(it.fw, it.rc, keys, vals, aux, prm, nk, cpos, fw_on);
        if (u < 0) { advance(it, r); continue; }
        u32 lo = (u32)u;
        u64 base = sshash ? words[prm[10] + lo] : (u64)aux[lo] + (u64)lo * (u64)(k - 1);
        int start = it.pos;
        if (start > last_pos) { emit(out, nh, hcap, start, lo, fw_on, cpos); last_pos = start; }
        i64 dir = fw_on ? 1 : -1;
        i64 dist = fw_on ? nk - 1 - cpos : cpos;
        i64 ccur = cpos;
        int lpos = start;
        i64 lcpos = cpos;
        bool lfw = fw_on, ended = false;
        while (!exhausted(it) && dist > 0) {
            i64 inc = (i64)advance(it, r);
            dist -= inc;
            if (dist < 0) break;
            ccur += dir * inc;
            if (exhausted(it)) break;
            u64 want = word_at(words, base + (u64)ccur, k);
            if (it.fw != want && it.rc != want) break;
            lfw = it.fw == want; lcpos += dir * inc; lpos = it.pos; ended = dist == 0;
        }
        if (lpos > last_pos) { emit(out, nh, hcap, lpos, lo, lfw, lcpos); last_pos = lpos; }
        if (ended || it.pos == start) advance(it, r);
    }
    counts[i] = nh;
}
"#;

/// The k-mer table on the GPU backend of this build.
pub struct KmerGpu {
    #[cfg(gpu)]
    dev: crate::gpu::KmerDevice,
}

impl KmerGpu {
    /// Whether this build has a GPU backend and it finds a device.
    pub fn available() -> bool {
        #[cfg(gpu)]
        {
            crate::gpu::available()
        }
        #[cfg(not(gpu))]
        {
            false
        }
    }

    /// Upload `table`; `Err` without a backend or device.
    pub fn new(table: &KmerTable) -> Result<Self, String> {
        #[cfg(gpu)]
        {
            let (keys, vals, aux, words, prm) = table.device();
            let dev = crate::gpu::KmerDevice::new(keys, vals, &aux, &words, &prm)?;
            Ok(Self { dev })
        }
        #[cfg(not(gpu))]
        {
            let _ = table;
            Err("this build has no GPU backend".into())
        }
    }

    /// The raw hits of the reads `seq[off[i]..off[i + 1]]`: read `i`'s first
    /// `min(counts[i], HITS_CAP)` hits are `hits[i * HITS_CAP..]`, and
    /// `counts[i] > HITS_CAP` means it had more. Safe to call from many threads.
    pub fn hits(
        &self,
        seq: &[u8],
        off: &[u32],
        hits: &mut Vec<Hit>,
        counts: &mut Vec<u32>,
    ) -> Result<(), String> {
        let n = off.len().saturating_sub(1);
        hits.clear();
        hits.resize(n * HITS_CAP, [0; 3]);
        counts.clear();
        counts.resize(n, 0);
        if n == 0 {
            return Ok(());
        }
        #[cfg(gpu)]
        {
            self.dev.hits(seq, off, hits.as_flattened_mut(), counts)
        }
        #[cfg(not(gpu))]
        {
            let _ = seq;
            Err("this build has no GPU backend".into())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Unitigs as plain strings, through the same `kmer(u, j)` callback the
    /// index uses.
    fn table(unitigs: &[&[u8]], k: usize) -> KmerTable {
        let mut kstart = vec![0u32];
        for u in unitigs {
            kstart.push(kstart.last().unwrap() + (u.len() - k + 1) as u32);
        }
        let kmer = |u: usize, j: u32| {
            let s = &unitigs[u][j as usize..j as usize + k];
            s.iter()
                .enumerate()
                .fold(0, |a, (i, &b)| a | enc(b) << (2 * i))
        };
        KmerTable::build(k, kstart, |_| true, kmer).unwrap()
    }

    fn hits(t: &KmerTable, read: &[u8]) -> Vec<Hit> {
        let mut v = Vec::new();
        t.hits(read, &mut v);
        v
    }

    fn rc_str(s: &[u8]) -> Vec<u8> {
        s.iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                x => x,
            })
            .collect()
    }

    /// Two 40-base pseudo-random unitigs (no 11-mer shared or repeated).
    fn unitigs() -> [Vec<u8>; 2] {
        let mut x = 0x9E37_79B9_u64;
        [0, 1].map(|_| {
            (0..40)
                .map(|_| {
                    x = x
                        .wrapping_mul(6364136223846793005)
                        .wrapping_add(1442695040888963407);
                    b"ACGT"[(x >> 62) as usize]
                })
                .collect()
        })
    }

    #[test]
    fn revcomp_matches_the_string_reverse_complement() {
        let enc_all = |s: &[u8]| {
            s.iter()
                .enumerate()
                .fold(0, |a, (i, &b)| a | enc(b) << (2 * i))
        };
        let [_, u1] = unitigs();
        for k in [5u32, 11, 31] {
            let s = &u1[..k as usize];
            assert_eq!(revcomp(enc_all(s), k), enc_all(&rc_str(s)));
        }
    }

    #[test]
    fn a_read_inside_a_unitig_is_one_walk() {
        let [u0, u1] = unitigs();
        let t = table(&[&u0, &u1], 11);
        // read = u1[3..20]: starts at k-mer 3 of unitig 1, walks to k-mer 9
        let h = hits(&t, &u1[3..20]);
        assert_eq!(h, [[0, 1 << 1 | 1, 3], [6, 1 << 1 | 1, 9]]);
        // its reverse complement walks the unitig backwards
        let h = hits(&t, &rc_str(&u1[3..20]));
        assert_eq!(h, [[0, 1 << 1, 9], [6, 1 << 1, 3]]);
    }

    #[test]
    fn substitution_splits_the_walk_and_n_is_jumped() {
        let [u0, u1] = unitigs();
        let t = table(&[&u0, &u1], 11);
        let mut read = u1[0..25].to_vec();
        read[12] = if read[12] == b'A' { b'C' } else { b'A' };
        let h = hits(&t, &read);
        // the walk from 0 stops at the first k-mer covering the substitution;
        // the next match is the first k-mer past it
        assert_eq!(h.first(), Some(&[0, 3, 0]));
        assert!(h.iter().any(|x| x[0] == 13 && x[2] == 13));
        // an N: windows covering it are skipped, the walk continues past it
        let mut read = u1[0..25].to_vec();
        read[12] = b'N';
        let h = hits(&t, &read);
        assert_eq!(h, [[0, 3, 0], [14, 3, 14]]);
    }

    #[test]
    fn short_reads_homopolymers_and_misses_give_no_hits() {
        let [u0, u1] = unitigs();
        let t = table(&[&u0, &u1], 11);
        assert!(hits(&t, &u1[..10]).is_empty());
        assert!(hits(&t, b"AAAAAAAAAAAA").is_empty());
        assert!(hits(&t, b"CCCCCCCCGGGGGGGG").is_empty());
        // lowercase is the same sequence
        assert_eq!(
            hits(&t, &u1[3..20].to_ascii_lowercase()),
            hits(&t, &u1[3..20])
        );
    }

    #[test]
    fn duplicate_kmers_and_bad_k_are_refused() {
        let kstart = vec![0, 1, 2];
        assert!(KmerTable::build(5, kstart.clone(), |_| true, |_, _| 0b1011).is_none());
        assert!(KmerTable::build(6, kstart.clone(), |_| true, |u, _| u as u64).is_none());
        assert!(KmerTable::build(33, kstart, |_| true, |u, _| u as u64).is_none());
    }

    /// The GPU kernel equals [`KmerTable::hits`] on random reads with errors,
    /// Ns and lowercase (skips without a device).
    #[cfg(gpu)]
    #[test]
    fn gpu_hits_match_the_cpu_reference() {
        use rand::{RngExt, SeedableRng};
        if !KmerGpu::available() {
            eprintln!("no GPU device; skipping");
            return;
        }
        let mut rng = rand_pcg::Pcg64Mcg::seed_from_u64(11);
        let bases = b"ACGT";
        let unitigs: Vec<Vec<u8>> = (0..300)
            .map(|_| {
                (0..rng.random_range(31..400))
                    .map(|_| bases[rng.random_range(0..4)])
                    .collect()
            })
            .collect();
        let refs: Vec<&[u8]> = unitigs.iter().map(|u| u.as_slice()).collect();
        let t = table(&refs, 31);
        let gpu = KmerGpu::new(&t).unwrap();
        let (mut seq, mut off) = (Vec::new(), vec![0u32]);
        let mut want = Vec::new();
        for i in 0..10_000 {
            let u = &unitigs[rng.random_range(0..unitigs.len())];
            let a = rng.random_range(0..u.len());
            let b = (a + rng.random_range(0..150)).min(u.len());
            let mut r = u[a..b].to_vec();
            if i % 3 == 0 {
                r = rc_str(&r);
            }
            for x in r.iter_mut() {
                match rng.random_range(0..200) {
                    0 => *x = bases[rng.random_range(0..4)],
                    1 => *x = b'N',
                    2 => x.make_ascii_lowercase(),
                    _ => {}
                }
            }
            want.push(hits(&t, &r));
            seq.extend_from_slice(&r);
            off.push(seq.len() as u32);
        }
        let (mut h, mut c) = (Vec::new(), Vec::new());
        gpu.hits(&seq, &off, &mut h, &mut c).unwrap();
        for (i, w) in want.iter().enumerate() {
            assert_eq!(c[i] as usize, w.len(), "read {i} count");
            let n = w.len().min(HITS_CAP);
            assert_eq!(&h[i * HITS_CAP..i * HITS_CAP + n], &w[..n], "read {i}");
        }
    }
}
