//! GPU k-mer lookup for pseudoalignment (`--sketch`), after kallisto's GPU
//! port (Melsted et al., bioRxiv 2026.03.04.709526): a flat open-addressing
//! table from canonical k-mer to its unitig position, probed for every k-mer of
//! every read by one GPU thread per read.
//!
//! The kernel reproduces piscem's hit collection exactly, with either skipping
//! strategy. `Strict` (`HitSearcher::walk_safely_until`) looks a k-mer up,
//! emits a raw hit, walks the read along that unitig while each next read
//! k-mer is the next unitig k-mer, and emits the walk's last matching k-mer as
//! a second hit. `Permissive` (`collect_permissive`, salmon's default) jumps
//! from a hit to the end of its unitig and checks only that k-mer, falling
//! back to a lookup there, one halfway, then the `Strict` walk.
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
    /// `fp << 32 | bucket` entries, [`EMPTY`] where free; empty when the
    /// device builds it (`cap` entries)
    table: Vec<u64>,
    cap: u32,
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
    // never all ones: the device build claims a slot with a CAS on it
    let fp = (h as u32).min(u32::MAX - 1);
    ((((h >> 32) * u64::from(cap)) >> 32) as u32, fp)
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
    /// `cpu_table`: build the minimizer → bucket table here, as [`Self::hits`]
    /// needs; without it [`KmerGpu::new`] builds it on the device (when
    /// [`KmerGpu::builds_table`]), and the result is only for the GPU.
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
        cpu_table: bool,
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
        if !cpu_table {
            let table = Vec::new();
            return Some(Self::sshash_table(
                k, m, magic, seq, soff, bits, table, cap, bstart, offs,
            ));
        }
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
        let table = table.into_iter().map(AtomicU64::into_inner).collect();
        unique.then(|| Self::sshash_table(k, m, magic, seq, soff, bits, table, cap, bstart, offs))
    }

    #[allow(clippy::too_many_arguments)]
    fn sshash_table(
        k: u32,
        m: u32,
        magic: u64,
        seq: Vec<u64>,
        soff: Vec<u64>,
        keep: Vec<u32>,
        table: Vec<u64>,
        cap: u32,
        bstart: Vec<u32>,
        offs: Vec<u32>,
    ) -> Self {
        Self {
            k,
            keys: Vec::new(),
            vals: Vec::new(),
            kstart: Vec::new(),
            seq,
            sshash: Some(Sshash {
                m,
                magic,
                soff,
                keep,
                table,
                cap,
                bstart,
                offs,
            }),
        }
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
        assert_eq!(
            s.table.len(),
            s.cap as usize,
            "the CPU lookup needs `cpu_table`"
        );
        let cap = s.cap;
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
    /// capacity, unitigs, hits per read, option B, m, magic (low, high),
    /// where `offs`, `keep` and `soff` start, and `Permissive` skipping.
    #[cfg_attr(not(gpu), allow(dead_code))]
    fn device(
        &self,
        permissive: bool,
    ) -> (&[u64], &[u32], Cow<'_, [u32]>, Cow<'_, [u64]>, [u32; 12]) {
        let hcap = HITS_CAP as u32;
        let perm = u32::from(permissive);
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
                perm,
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
            s.cap,
            s.soff.len() as u32 - 1,
            hcap,
            1,
            s.m,
            s.magic as u32,
            (s.magic >> 32) as u32,
            s.bstart.len() as u32,
            (s.bstart.len() + s.offs.len()) as u32,
            self.seq.len() as u32,
            perm,
        ];
        // `vals` is unused here; one word, as CUDA refuses an empty upload
        (&s.table, &[0], Cow::Owned(aux), Cow::Owned(words), prm)
    }

    /// The raw hits of `read`, appended to `out`, as piscem's `Permissive`
    /// (`permissive`) or `Strict` skipping finds them: the CPU reference of the
    /// GPU kernel ([`KMER_BODY`]), which computes the same in one loop.
    pub fn hits(&self, read: &[u8], permissive: bool, out: &mut Vec<Hit>) {
        let k = self.k as i32;
        let len = read.len() as i32;
        let mut it = Iter::new(read, k);
        let mut last = -1;
        if !permissive {
            self.walk(&mut it, len - k, &mut last, out);
            return;
        }
        // piscem's `collect_permissive`: from a hit, jump to the end of its
        // unitig (or of the read) and check that k-mer on the unitig text;
        // failing that, a lookup there, then one halfway, then a safe walk
        while !it.exhausted() {
            if it.homopolymer() {
                it.advance();
                continue;
            }
            let Some((u, nk, cpos, fw_on)) = self.locate(it.fw, it.rc) else {
                it.advance();
                continue;
            };
            let pos = it.pos;
            if pos > last {
                push(out, &mut last, pos, u, fw_on, cpos);
            }
            let dir: i64 = if fw_on { 1 } else { -1 };
            let skip =
                i64::from(len - k - pos).min(if fw_on { nk - 1 - cpos } else { cpos }) as i32;
            if skip <= 1 {
                it.advance();
                continue;
            }
            let backup = it.clone();
            let nd = it.advance();
            if it.exhausted() {
                continue;
            }
            if nd < skip {
                // the next k-mer leaves the unitig: look it up instead
                if !self.direct(&it, u, fw_on, cpos + dir * i64::from(nd)) {
                    continue;
                }
                it = backup.clone();
            }
            let moved = it.advance_by(skip);
            if it.exhausted() {
                it = backup;
                it.advance();
                self.walk(&mut it, pos + skip, &mut last, out);
                continue;
            }
            let target = it.clone();
            let c = cpos + dir * i64::from(skip);
            if moved == skip && self.direct(&it, u, fw_on, c) {
                push(out, &mut last, it.pos, u, fw_on, c);
                it.advance();
                continue;
            }
            let alt = self.locate(it.fw, it.rc);
            if let Some((u2, _, c2, f2)) = alt {
                if u2 == u && f2 == fw_on && if fw_on { c2 > cpos } else { c2 < cpos } {
                    push(out, &mut last, it.pos, u2, f2, c2);
                    it.advance();
                    continue;
                }
            }
            let mut mid_ok = false;
            if skip > 4 {
                let mut mid = backup.clone();
                mid.advance_by(skip / 2);
                if !mid.exhausted() {
                    if let Some((um, _, cm, fm)) = self.locate(mid.fw, mid.rc) {
                        if um == u {
                            push(out, &mut last, mid.pos, um, fm, cm);
                            if let Some((u2, _, c2, f2)) = alt {
                                push(out, &mut last, target.pos, u2, f2, c2);
                            }
                            mid_ok = true;
                        } else if let Some((u2, _, c2, f2)) = alt.filter(|a| a.0 == um) {
                            push(out, &mut last, target.pos, u2, f2, c2);
                            mid_ok = true;
                        }
                    }
                }
            }
            if mid_ok {
                it = target;
                it.advance();
            } else {
                it = backup;
                it.advance();
                it.advance();
                self.walk(&mut it, target.pos, &mut last, out);
            }
        }
    }

    /// Whether the read k-mer at `it` is unitig `u`'s k-mer `c`, in the
    /// orientation `fw_on` (piscem's `check_direct_match`).
    fn direct(&self, it: &Iter, u: usize, fw_on: bool, c: i64) -> bool {
        let want = word_at(&self.seq, self.base(u) + c as u64, self.k);
        let fw = it.fw == want;
        (fw || it.rc == want) && fw == fw_on
    }

    /// piscem's `walk_safely_until`: every k-mer from `it` up to read position
    /// `end`, a hit extended along its unitig for as long as the read matches it.
    fn walk(&self, it: &mut Iter, end: i32, last: &mut i32, out: &mut Vec<Hit>) {
        while !it.exhausted() && it.pos <= end {
            if it.homopolymer() {
                it.advance();
                continue;
            }
            let Some((u, nk, cpos, fw_on)) = self.locate(it.fw, it.rc) else {
                it.advance();
                continue;
            };
            let base = self.base(u);
            let start = it.pos;
            if start > *last {
                push(out, last, start, u, fw_on, cpos);
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
            if lpos > *last {
                push(out, last, lpos, u, lfw, lcpos);
            }
            if ended || it.pos == start {
                it.advance();
            }
        }
    }
}

fn push(out: &mut Vec<Hit>, last: &mut i32, pos: i32, u: usize, fw: bool, cpos: i64) {
    out.push([pos as u32, (u as u32) << 1 | u32::from(fw), cpos as u32]);
    *last = pos;
}

/// piscem's `ReadKmerIter`: the read's k-mer windows free of non-ACGT bytes.
#[derive(Clone)]
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

    fn homopolymer(&self) -> bool {
        let w = &self.read[self.pos as usize..(self.pos + self.k) as usize];
        w.iter().all(|&b| b == w[0])
    }

    /// Move `n` windows right (past any non-ACGT byte); the distance moved.
    fn advance_by(&mut self, n: i32) -> i32 {
        if n <= 0 || self.exhausted() {
            return 0;
        }
        let old = self.pos;
        if n <= self.k {
            for _ in 0..n {
                self.advance();
                if self.exhausted() {
                    break;
                }
            }
        } else {
            // a long jump rebuilds the window from scratch
            self.pos += n;
            let len = self.read.len() as i32;
            if self.pos + self.k > len {
                return len - old;
            }
            if let Some(i) = (self.pos..self.pos + self.k).rfind(|&i| !valid(self.read[i as usize]))
            {
                self.last_inv = i;
            }
            if self.last_inv >= self.pos {
                self.pos = self.last_inv + 1;
                self.find_next_valid();
            }
            if !self.exhausted() {
                self.build();
            }
        }
        self.pos - old
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
FN void emit(G u32* out, PRIV u32& nh, PRIV int& last, u32 hcap, int pos, u32 u, bool fw, i64 cpos) {
    if (nh < hcap) {
        out[3 * nh] = (u32)pos;
        out[3 * nh + 1] = (u << 1) | (fw ? 1u : 0u);
        out[3 * nh + 2] = (u32)cpos;
    }
    nh += 1;
    last = pos;
}
FN bool homopolymer(PRIV const Iter& it, G const u8* r) {
    for (int j = 1; j < it.k; ++j) {
        if (r[it.pos + j] != r[it.pos]) return false;
    }
    return true;
}
FN int advance_by(PRIV Iter& it, G const u8* r, int n) {
    if (n <= 0 || exhausted(it)) return 0;
    int old = it.pos;
    if (n <= it.k) {
        for (int j = 0; j < n; ++j) {
            advance(it, r);
            if (exhausted(it)) break;
        }
    } else {
        it.pos += n;
        if (it.pos + it.k > it.len) return it.len - old;
        for (int j = it.pos; j < it.pos + it.k; ++j) {
            if (!valid_base(r[j])) it.last_inv = j;
        }
        if (it.last_inv >= it.pos) { it.pos = it.last_inv + 1; find_next_valid(it, r); }
        if (!exhausted(it)) build_words(it, r);
    }
    return it.pos - old;
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
    u32 s = (u32)(((h >> 32) * (u64)cap) >> 32), fp = (u32)h < 0xFFFFFFFFu ? (u32)h : 0xFFFFFFFEu;
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
// the lookup's buffers, as parameters and as arguments
#define CTX G const u64* keys, G const u32* vals, G const u32* aux, G const u64* words, G const u32* prm
#define CTXA keys, vals, aux, words, prm
FN i64 locate(u64 fw, u64 rc, CTX, PRIV i64& nk, PRIV i64& cpos, PRIV bool& fw_on) {
    return prm[4] != 0 ? locate_b(fw, rc, keys, aux, words, prm, nk, cpos, fw_on)
                       : locate_a(fw, rc, keys, vals, aux, prm, nk, cpos, fw_on);
}
// unitig `u`'s first base in `words`
FN u64 unitig_base(u32 u, CTX) {
    return prm[4] != 0 ? words[prm[10] + u] : (u64)aux[u] + (u64)u * (u64)(prm[0] - 1);
}
// `KmerTable::direct`
FN bool direct(PRIV const Iter& it, u64 base, bool fw_on, i64 c, CTX) {
    u64 want = word_at(words, base + (u64)c, (int)prm[0]);
    bool fw = it.fw == want;
    return (fw || it.rc == want) && fw == fw_on;
}
// prm: see `KmerTable::device`. `KmerTable::hits` as one loop making one
// lookup per turn, so the lookup and the walk are each inlined once (the
// kernel's size sets its occupancy). `walk_end`: the `Strict` walk runs up to
// that read position; `phase`: what the lookup is for (0: a walk or a new
// anchor, 1: the k-mer a skip landed on, 2: the one halfway).
FN void map_mate(u32 i, G const u8* seq, G const u32* off, CTX, G u32* hits, G u32* counts) {
    G const u8* r = seq + off[i];
    int k = (int)prm[0];
    u32 hcap = prm[3];
    G u32* out = hits + (u64)i * hcap * 3;
    u32 nh = 0;
    int last = -1;
    Iter it;
    it.pos = 0; it.last_inv = -1; it.len = (int)(off[i + 1] - off[i]); it.k = k; it.fw = 0; it.rc = 0;
    find_next_valid(it, r);
    if (!exhausted(it)) build_words(it, r);
    int walk_end = prm[11] == 0 ? it.len : -1;
    int phase = 0;
    // the anchor of a skip, and the lookup at its target
    Iter backup = it, mid = it;
    u32 au = 0;
    bool afw = false;
    i64 acpos = 0, u2 = -1, c2 = 0;
    bool f2 = false;
    int skip = 0;
    for (;;) {
        if (phase == 0) {
            if (exhausted(it)) break;
            if (it.pos > walk_end) walk_end = -1;
            if (homopolymer(it, r)) { advance(it, r); continue; }
        }
        bool fw_on = false;
        i64 nk = 0, cpos = 0;
        u64 qf = phase == 2 ? mid.fw : it.fw, qr = phase == 2 ? mid.rc : it.rc;
        i64 u = locate(qf, qr, CTXA, nk, cpos, fw_on);
        if (phase == 1) {
            u2 = u; c2 = cpos; f2 = fw_on;
            if (u == (i64)au && f2 == afw && (afw ? c2 > acpos : c2 < acpos)) {
                emit(out, nh, last, hcap, it.pos, au, f2, c2);
                advance(it, r);
                phase = 0;
                continue;
            }
            if (skip > 4) {
                mid = backup;
                advance_by(mid, r, skip / 2);
                if (!exhausted(mid)) { phase = 2; continue; }
            }
        } else if (phase == 2) {
            bool ok = false;
            if (u >= 0 && u == (i64)au) {
                emit(out, nh, last, hcap, mid.pos, au, fw_on, cpos);
                if (u2 >= 0) emit(out, nh, last, hcap, it.pos, (u32)u2, f2, c2);
                ok = true;
            } else if (u >= 0 && u == u2) {
                emit(out, nh, last, hcap, it.pos, (u32)u2, f2, c2);
                ok = true;
            }
            if (ok) { advance(it, r); phase = 0; continue; }
        }
        if (phase != 0) {
            // neither the target nor halfway confirms the skip: walk to it
            int target = it.pos;
            it = backup;
            advance(it, r);
            advance(it, r);
            walk_end = target;
            phase = 0;
            continue;
        }
        if (u < 0) { advance(it, r); continue; }
        u32 lo = (u32)u;
        u64 base = unitig_base(lo, CTXA);
        int pos = it.pos;
        if (pos > last) emit(out, nh, last, hcap, pos, lo, fw_on, cpos);
        i64 dir = fw_on ? 1 : -1;
        if (walk_end >= 0) {
            // `KmerTable::walk`
            i64 dist = fw_on ? nk - 1 - cpos : cpos;
            i64 ccur = cpos;
            int lpos = pos;
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
            if (lpos > last) emit(out, nh, last, hcap, lpos, lo, lfw, lcpos);
            if (ended || it.pos == pos) advance(it, r);
            continue;
        }
        // `Permissive`: skip to the unitig's or the read's end
        i64 dc = fw_on ? nk - 1 - cpos : cpos, dr = (i64)(it.len - k - pos);
        skip = (int)(dr < dc ? dr : dc);
        if (skip <= 1) { advance(it, r); continue; }
        backup = it;
        int nd = advance(it, r);
        if (exhausted(it)) continue;
        if (nd < skip) {
            if (!direct(it, base, fw_on, cpos + dir * (i64)nd, CTXA)) continue;
            it = backup;
        }
        int moved = advance_by(it, r, skip);
        if (exhausted(it)) {
            it = backup;
            advance(it, r);
            walk_end = pos + skip;
            continue;
        }
        i64 c = cpos + dir * (i64)skip;
        if (moved == skip && direct(it, base, fw_on, c, CTXA)) {
            emit(out, nh, last, hcap, it.pos, lo, fw_on, c);
            advance(it, r);
            continue;
        }
        au = lo; afw = fw_on; acpos = cpos;
        phase = 1;
    }
    counts[i] = nh;
}
#undef CTX
#undef CTXA
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

    /// Whether this backend builds the sshash lookup's minimizer table on the
    /// device (`KmerTable::build_sshash` with `cpu_table` off): Metal.
    pub fn builds_table() -> bool {
        cfg!(gpu = "metal")
    }

    /// Upload `table`, for piscem's `Permissive` (`permissive`) or `Strict`
    /// skipping; `Err` without a backend or device.
    pub fn new(table: &KmerTable, permissive: bool) -> Result<Self, String> {
        #[cfg(gpu)]
        {
            let (keys, vals, aux, words, prm) = table.device(permissive);
            let dev = crate::gpu::KmerDevice::new(keys, vals, &aux, &words, &prm)?;
            Ok(Self { dev })
        }
        #[cfg(not(gpu))]
        {
            let _ = (table, permissive);
            Err("this build has no GPU backend".into())
        }
    }

    /// Run the kernel until the device is up to speed: its first calls take
    /// up to ten times longer (clocks ramping up, buffers first touched).
    /// The reads are cut from `table`'s unitig text.
    pub fn warm_up(&self, table: &KmerTable) {
        const LEN: u64 = 100;
        let bases = (table.seq.len() as u64).saturating_sub(2) * 32;
        if bases <= LEN {
            return;
        }
        let (mut seq, mut off) = (Vec::new(), vec![0u32]);
        let mut x = 0x9E37_79B9_7F4A_7C15u64;
        for _ in 0..8192 {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let p = (x >> 11) % (bases - LEN);
            seq.extend(
                (p..p + LEN).map(|b| {
                    b"ACTG"[(table.seq[(b / 32) as usize] >> (2 * (b % 32)) & 3) as usize]
                }),
            );
            off.push(seq.len() as u32);
        }
        let (mut h, mut c) = (Vec::new(), Vec::new());
        let mut best = f64::MAX;
        for i in 0..40 {
            let t = std::time::Instant::now();
            if self.hits(&seq, &off, &mut h, &mut c).is_err() {
                return;
            }
            let dt = t.elapsed().as_secs_f64();
            best = best.min(dt);
            if i >= 8 && dt < 1.25 * best {
                return;
            }
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

    /// `Strict` hits.
    fn hits(t: &KmerTable, read: &[u8]) -> Vec<Hit> {
        let mut v = Vec::new();
        t.hits(read, false, &mut v);
        v
    }

    fn permissive(t: &KmerTable, read: &[u8]) -> Vec<Hit> {
        let mut v = Vec::new();
        t.hits(read, true, &mut v);
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
    fn permissive_jumps_to_the_read_end_over_a_substitution() {
        let [u0, u1] = unitigs();
        let t = table(&[&u0, &u1], 11);
        assert_eq!(permissive(&t, &u1[3..20]), hits(&t, &u1[3..20]));
        // the k-mer at the read's end is checked on the unitig, the
        // substitution between is never seen (`Strict` splits the walk there)
        let mut read = u1[0..25].to_vec();
        read[12] = if read[12] == b'A' { b'C' } else { b'A' };
        assert_eq!(permissive(&t, &read), [[0, 3, 0], [14, 3, 14]]);
        assert_ne!(hits(&t, &read), [[0, 3, 0], [14, 3, 14]]);
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
        let (mut seq, mut off, mut reads) = (Vec::new(), vec![0u32], Vec::new());
        for i in 0..10_000 {
            let u = &unitigs[rng.random_range(0..unitigs.len())];
            let a = rng.random_range(0..u.len());
            let b = (a + rng.random_range(0..150)).min(u.len());
            let mut r = u[a..b].to_vec();
            if i % 3 == 0 {
                r = rc_str(&r);
            }
            // an occasional jump to another unitig, so skips meet mismatches
            if i % 4 == 1 {
                let v = &unitigs[rng.random_range(0..unitigs.len())];
                r.extend_from_slice(&v[..v.len().min(60)]);
            }
            for x in r.iter_mut() {
                match rng.random_range(0..200) {
                    0 => *x = bases[rng.random_range(0..4)],
                    1 => *x = b'N',
                    2 => x.make_ascii_lowercase(),
                    _ => {}
                }
            }
            seq.extend_from_slice(&r);
            off.push(seq.len() as u32);
            reads.push(r);
        }
        for perm in [false, true] {
            let gpu = KmerGpu::new(&t, perm).unwrap();
            let (mut h, mut c) = (Vec::new(), Vec::new());
            gpu.hits(&seq, &off, &mut h, &mut c).unwrap();
            for (i, r) in reads.iter().enumerate() {
                let mut w = Vec::new();
                t.hits(r, perm, &mut w);
                assert_eq!(c[i] as usize, w.len(), "permissive {perm} read {i} count");
                let n = w.len().min(HITS_CAP);
                assert_eq!(
                    &h[i * HITS_CAP..i * HITS_CAP + n],
                    &w[..n],
                    "permissive {perm} read {i}"
                );
            }
        }
    }
}
