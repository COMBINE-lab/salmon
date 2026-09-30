# Unspliced targets on public total RNA-seq

Benchmark of `salmon index --unspliced` on public rRNA-depleted RNA-seq,
run with `scripts/bench_unspliced.sh` and summarized with
`scripts/bench_unspliced.py` (outputs in this directory).

## Data and setup

- **Reads**: GSE100127 (ABRF rRNA-depletion study, Herbert et al. 2018,
  [doi:10.1186/s12864-018-4585-1](https://doi.org/10.1186/s12864-018-4585-1)).
  Universal Human Reference RNA, Ribo-Zero Gold, reverse-stranded, 2x76,
  technical duplicate libraries of the same RNA from two sites:
  SRR5693352 / SRR5693353 (site 2), SRR5693399 / SRR5693400 (site 4);
  10.3 to 16.0 M read pairs each.
- **References**: GENCODE v50 transcripts (`--gencode`), GRCh38 primary
  assembly as decoys, `gencode.v50.primary_assembly.annotation.gtf.gz`.
- **Indices** (k = 31, 16 threads): `gentrome`; `intron` (`--unspliced intron
  --readLength 76`, flank 75); `premrna` (`--unspliced premrna`). Each
  unspliced index was built in both layouts: `projection` (default with
  genome decoys) and `sequence` (`--unsplicedLayout sequence`).
- **Quant**: `-l A --gcBias --seqBias -p 16`, deterministic default path.
- **Machine**: Apple M-series, 16 cores, 128 GB, shared with unrelated jobs
  during the runs (load average 50 to 100). Absolute times are noisy; the
  ratios below compare runs made in the same series, where all modes shared
  the same conditions.

## Where the fragments go (projection layout)

| mode | mapping rate | to decoy | on unspliced targets | retained_intron share of transcript reads |
|---|---|---|---|---|
| gentrome | 0.641 to 0.696 | 0.239 to 0.305 | 0.010 | 0.043 to 0.047 |
| intron | 0.912 to 0.916 | 0.015 to 0.016 | 0.255 to 0.312 | 0.029 to 0.031 |
| premrna | 0.913 to 0.918 | 0.013 to 0.015 | 0.263 to 0.319 | 0.030 to 0.032 |

(The ~1% "on unspliced targets" in the gentrome rows is transcripts of the
FASTA that the primary-assembly GTF does not list.) With intron targets, 91% of fragments are
assigned instead of 64 to 70%; what remains on the decoy is intergenic or
aligns better elsewhere. The share of transcript reads on `retained_intron`
isoforms drops by about a third.

`quant.genes.usa.tsv` / `meta_info.json`, `intron` mode: 53 to 60% of assigned
fragments spliced, 24 to 30% unspliced, 16 to 17% ambiguous; `premrna` moves
exonic fragments of intron-bearing genes to ambiguous (37 to 45% spliced,
24 to 30% unspliced, 31 to 33% ambiguous). Decoy fragments are 1.4 to 1.7% of
assigned + decoy.

## Reproducibility of gene avgTxLength between technical replicates

Gene-level average transcript length (TPM-weighted mean effective length of
the spliced isoforms, as tximport computes it), for the 5,024 genes with more
than 250 spliced reads in every library under every mode (projection layout):

| mode | SD log2 ratio, SRR5693352 / SRR5693353 | SD log2 ratio, SRR5693399 / SRR5693400 | r of the site difference, replicate a vs b |
|---|---|---|---|
| gentrome | 0.139 | 0.195 | 0.202 |
| intron | 0.168 | 0.224 | 0.191 |
| premrna | 0.169 | 0.226 | 0.188 |

The sequence layout gives 0.158 / 0.208 (intron) and 0.159 / 0.211 (premrna)
on its own common gene set (`report-sequence.md`).

**No improvement here, rather a small loss.** This data set cannot show the
effect the option targets: within each replicate pair the pre-mRNA content is
nearly the same (30.5% vs 28.8% and 23.9% vs 26.4% of fragments to the decoy
under the gentrome), so the confound between pre-mRNA content and isoform
proportions is small, while the extra targets add EM ambiguity (next
section). A cohort whose libraries differ in pre-mRNA content is the relevant
test; with n = 4 the correlation between retained-intron share and decoy share
(1.00 gentrome, 0.80 intron or premrna) is not interpretable either.

## Projection vs sequence layout

On the simulation in the test suite the two layouts agree to 0.0001% of the
fragment mass. On SRR5693352 (`intron`), gene totals agree to 2.6% of the
mass, but a few pairs of targets the EM cannot tell apart land in different
corners of a flat likelihood: MALAT1 is 64.6 k spliced / 0.6 k unspliced in
the sequence layout and 5.9 k / 59.3 k in the projection layout (its intron
target covers what its single-exon isoform covers), and 81.7 k fragments swap
between RN7SK (ENSG00000202198) and ENSG00000283293. The equivalence classes
themselves differ slightly (835,320 vs 831,384), from how the mapper caps
repeated hits when intronic k-mers occur twice (sequence) or once
(projection).

## Cost

Quant wall time, summed over the four libraries, each series run under its
own load:

| series | gentrome | intron | premrna |
|---|---|---|---|
| sequence layout | 983 s | 1,299 s (+32%) | 1,516 s (+54%) |
| projection layout | 1,259 s | 1,228 s (-2%) | 1,195 s (-5%) |

Peak quant memory: 13.4 GB (gentrome), 18.4 to 18.5 GB (projection), 18.9 to
19.0 GB (sequence); the reference store holds the unspliced sequences, which
effective lengths and bias correction read.

Index: gentrome 9.8 GiB on disk, 422 s, 24.7 GB peak; intron 12.1 GiB
(projection) or 14.2 GiB (sequence), premrna 12.2 or 14.3 GiB; builds took
550 to 640 s at 31.5 to 38.4 GB peak (projection builds overlapped another
build on the same machine).

Profiling one library along the way (each comparison run back to back): with
the sequence layout and exact bias sweeps, intron targets cost more than twice
the gentrome's quant time; sampling unspliced targets at every 16th start in
the bias sweeps brought that to +41% (277 s vs 197 s), the rest being mapping,
which the projection layout removes.
