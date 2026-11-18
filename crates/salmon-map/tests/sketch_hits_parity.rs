//! The GPU k-mer lookup path (`salmon_infer::kmers` raw hits + the replay in
//! `salmon_map::sketch_hits`) must map every read exactly as piscem's CPU path
//! does, with `SkippingStrategy::Permissive` and with `Strict`. Runs without a GPU: the hits come
//! from `KmerTable::hits`, the CPU reference the GPU kernel is tested against.

use piscem_rs::mapping::hit_searcher::{HitSearcher, SkippingStrategy};
use piscem_rs::mapping::streaming_query::PiscemStreamingQuery;
use salmon_index::{build, IndexBuildOptions, SalmonIndex};
use salmon_infer::kmers::{Hit, KmerGpu, KmerTable, HITS_CAP};
use salmon_map::sketch_hits::{
    kept, map_read_pair_hits_into, map_single_read_hits_into, sshash_parts, unitig_kmers,
};
use salmon_map::{
    map_read_pair_sketch_into, map_single_read_sketch_into, RefProvider, SketchScratch,
};
use std::fmt::Write as _;

struct IdxRefs<'a>(&'a SalmonIndex);
impl RefProvider for IdxRefs<'_> {
    fn num_refs(&self) -> usize {
        self.0.num_refs()
    }
    fn ref_view(&self, tid: u32) -> salmon_core::RefView<'_> {
        salmon_core::RefView::Bytes(self.0.ref_seq(tid))
    }
}

struct Lcg(u64);
impl Lcg {
    fn next(&mut self, n: usize) -> usize {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 33) % n as u64) as usize
    }
    fn seq(&mut self, n: usize) -> Vec<u8> {
        (0..n).map(|_| b"ACGT"[self.next(4)]).collect()
    }
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

/// Transcripts sharing exons (so unitigs split and k-mers multimap), one
/// repeated element in many of them (so some hits are too repetitive), and
/// a palindromic stretch (so walks flip orientation).
fn transcripts(rng: &mut Lcg) -> Vec<Vec<u8>> {
    let exons: Vec<Vec<u8>> = (0..12)
        .map(|_| {
            let n = 60 + rng.next(200);
            rng.seq(n)
        })
        .collect();
    let repeat = rng.seq(45);
    let half = rng.seq(40);
    let pal = [half.clone(), revcomp(&half)].concat();
    (0..40)
        .map(|t| {
            let mut s = Vec::new();
            for (e, exon) in exons.iter().enumerate() {
                if (t * 7 + e * 3) % 5 < 3 {
                    s.extend_from_slice(exon);
                }
                if e == 4 && t % 3 == 0 {
                    s.extend_from_slice(&repeat);
                }
                if e == 8 && t % 11 == 0 {
                    s.extend_from_slice(&pal);
                }
            }
            s.extend(rng.seq(30));
            s
        })
        .collect()
}

/// A read of `len` from `t` at `at`, with substitutions, an occasional N and
/// lowercase base.
fn read(rng: &mut Lcg, t: &[u8], at: usize, len: usize) -> Vec<u8> {
    let mut r = t[at..(at + len).min(t.len())].to_vec();
    for b in r.iter_mut() {
        match rng.next(200) {
            0 | 1 => *b = b"ACGT"[rng.next(4)],
            2 => *b = b'N',
            3 => *b = b.to_ascii_lowercase(),
            _ => {}
        }
    }
    r
}

#[test]
fn hits_path_maps_like_piscem_strict() {
    let mut rng = Lcg(0x5A1_C0FFEE);
    let txs = transcripts(&mut rng);
    let tmp = tempfile::tempdir().unwrap();
    let fasta = tmp.path().join("tx.fa");
    let mut fa = String::new();
    for (i, t) in txs.iter().enumerate() {
        writeln!(fa, ">t{i}\n{}", String::from_utf8_lossy(t)).unwrap();
    }
    std::fs::write(&fasta, fa).unwrap();
    let out = tmp.path().join("idx");
    let mut opts = IndexBuildOptions::new(vec![fasta], out.clone());
    opts.threads = 1;
    build(&opts).expect("build index");
    let idx = SalmonIndex::load(&out).expect("load index");
    let index = idx.inner();
    let refs = IdxRefs(&idx);
    assert!(index.has_ec_table(), "the ambiguous-class path must run");

    let (kstart, keep, kmer) = unitig_kmers(index);
    let table = KmerTable::build(index.k(), kstart, keep, kmer).expect("table");
    // option B: piscem's own sshash dictionary, minimizer → bucket → positions
    let p = sshash_parts(index).expect("u32 offsets");
    let sshash = KmerTable::build_sshash(
        index.k(),
        p.m,
        p.magic,
        p.spss,
        p.soff,
        kept(index),
        p.bstart,
        p.offs,
        true,
    )
    .expect("sshash table");
    let mut hs = HitSearcher::new(index);
    let (mut cpu_s, mut hit_s) = (SketchScratch::new(index.k()), SketchScratch::new(index.k()));
    let (mut want, mut got) = (Vec::new(), Vec::new());
    let (mut h1, mut h2): (Vec<Hit>, Vec<Hit>) = (Vec::new(), Vec::new());
    let (mut mapped, mut differ) = (0, 0);
    // every read, packed as the GPU takes them
    let (mut seq, mut off) = (Vec::new(), vec![0u32]);

    for n in 0..3000 {
        let t = &txs[rng.next(txs.len())];
        let flen = (150 + rng.next(250)).min(t.len());
        let at = rng.next(t.len() - flen + 1);
        let rl = 20 + rng.next(130);
        let mut r1 = read(&mut rng, t, at, rl.min(flen));
        let mut r2 = revcomp(&read(&mut rng, t, at + flen - rl.min(flen), rl.min(flen)));
        if n % 2 == 1 {
            std::mem::swap(&mut r1, &mut r2);
        }
        if n % 17 == 0 {
            r2 = rng.seq(rl);
        }
        // small caps too, so the repetitive-hit, recovery and class paths run
        let (max_hit, max_read) = if n % 3 == 0 { (3, 4) } else { (256, 2500) };

        // the oracle's raw hits are piscem's, with either skipping strategy
        for r in [&r1, &r2] {
            seq.extend_from_slice(r);
            off.push(seq.len() as u32);
            let mut strict = Vec::new();
            for strat in [SkippingStrategy::Strict, SkippingStrategy::Permissive] {
                let perm = strat == SkippingStrategy::Permissive;
                sshash_lib::dispatch_on_k!(index.k(), K => {
                    let mut q = PiscemStreamingQuery::<K>::new(index.dict());
                    hs.get_raw_hits_sketch::<K>(r, &mut q, strat, true);
                });
                let piscem: Vec<Hit> = hs
                    .left_hits()
                    .iter()
                    .map(|(p, h)| {
                        [
                            *p as u32,
                            h.contig_id() << 1 | h.hit_fw_on_contig() as u32,
                            h.contig_pos(),
                        ]
                    })
                    .collect();
                for (name, t) in [("table", &table), ("sshash", &sshash)] {
                    h1.clear();
                    t.hits(r, perm, &mut h1);
                    assert_eq!(
                        h1,
                        piscem,
                        "{name} {strat:?} raw hits of {}",
                        String::from_utf8_lossy(r)
                    );
                }
                if perm {
                    differ += (piscem != strict) as usize;
                } else {
                    strict = piscem;
                }
            }
        }

        for strat in [SkippingStrategy::Strict, SkippingStrategy::Permissive] {
            let perm = strat == SkippingStrategy::Permissive;
            h1.clear();
            h2.clear();
            table.hits(&r1, perm, &mut h1);
            table.hits(&r2, perm, &mut h2);
            map_read_pair_sketch_into(
                &mut want, index, &mut hs, &mut cpu_s, &r1, &r2, false, false, strat, max_hit,
                max_read, &refs, false,
            );
            assert!(map_read_pair_hits_into(
                &mut got,
                index,
                &mut hit_s,
                (r1.len(), &h1),
                (r2.len(), &h2),
                false,
                false,
                max_hit,
                max_read,
                &refs,
                false,
            ));
            assert_eq!(got, want, "{strat:?} pair {n}");
            mapped += (perm && !want.is_empty()) as usize;

            map_single_read_sketch_into(
                &mut want, index, &mut hs, &mut cpu_s, &r1, strat, max_hit, max_read,
            );
            assert!(map_single_read_hits_into(
                &mut got,
                index,
                &mut hit_s,
                r1.len(),
                &h1,
                max_hit,
                max_read
            ));
            assert_eq!(got, want, "{strat:?} single {n}");
        }
    }
    assert!(mapped > 1500, "only {mapped} pairs mapped");
    // the skips of `Permissive` are exercised, not only its `Strict` fallbacks
    assert!(
        differ > 500,
        "Permissive differs from Strict on only {differ} reads"
    );
    eprintln!(
        "Permissive raw hits differ from Strict on {differ} of {} reads",
        off.len() - 1
    );

    // with a device, the kernel matches the CPU reference in both modes
    if !KmerGpu::available() {
        eprintln!("no GPU device; skipping the kernel check");
        return;
    }
    // and so does the minimizer table built on the device, where it is
    let on_device = KmerGpu::builds_table().then(|| {
        let p = sshash_parts(index).expect("u32 offsets");
        let (k, spss) = (index.k(), p.spss);
        KmerTable::build_sshash(
            k,
            p.m,
            p.magic,
            spss,
            p.soff,
            kept(index),
            p.bstart,
            p.offs,
            false,
        )
        .expect("sshash parts")
    });
    let mut runs = vec![("table", &table, &table), ("sshash", &sshash, &sshash)];
    runs.extend(
        on_device
            .iter()
            .map(|d| ("sshash, device-built table", d, &sshash)),
    );
    for (name, t, cpu) in runs {
        for perm in [false, true] {
            let gpu = KmerGpu::new(t, perm).expect("upload");
            gpu.warm_up(t);
            let (mut h, mut c) = (Vec::new(), Vec::new());
            gpu.hits(&seq, &off, &mut h, &mut c).expect("kernel");
            for i in 0..off.len() - 1 {
                h1.clear();
                cpu.hits(&seq[off[i] as usize..off[i + 1] as usize], perm, &mut h1);
                assert_eq!(
                    c[i] as usize,
                    h1.len(),
                    "{name} permissive {perm} read {i} count"
                );
                let n = h1.len().min(HITS_CAP);
                assert_eq!(
                    &h[i * HITS_CAP..i * HITS_CAP + n],
                    &h1[..n],
                    "{name} permissive {perm} read {i}"
                );
            }
        }
    }
}
