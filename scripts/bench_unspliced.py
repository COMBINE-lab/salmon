#!/usr/bin/env python3
"""Summarize a `scripts/bench_unspliced.sh` run.

    scripts/bench_unspliced.py OUT_DIR GTF [--pairs A,B C,D ...] [--min-reads 250]

Per index mode and library: mapping rate, share of fragments dropped as decoy,
share on unspliced targets, and the share of *spliced* (transcript) reads that
land on `retained_intron` isoforms. Then, per mode, the reproducibility of
isoform proportions between technical replicate libraries, measured the way
tximport uses them: the gene-level average transcript length (TPM-weighted
mean of the isoforms' effective lengths), for genes whose spliced reads exceed
--min-reads in every library of the comparison. Timings come from times.tsv.

Writes OUT_DIR/summary.tsv (per mode and library) and prints Markdown tables.
Needs numpy, pandas and scipy.
"""

import argparse
import gzip
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


def gtf_transcripts(path):
    """transcript_id -> (gene_id, transcript_type), from the GTF's transcript lines."""
    rows = {}
    attr = re.compile(r'(\w+) "([^"]*)"')
    with gzip.open(path, "rt") if str(path).endswith(".gz") else open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "transcript":
                continue
            a = dict(attr.findall(cols[8]))
            rows[a["transcript_id"]] = (a["gene_id"], a.get("transcript_type", ""))
    return pd.DataFrame.from_dict(rows, orient="index", columns=["gene", "type"])


def load(out, mode, run, tx):
    q = out / "quant" / mode / run
    sf = pd.read_csv(q / "quant.sf", sep="\t", index_col=0)
    meta = json.loads((q / "aux_info" / "meta_info.json").read_text())
    spliced = sf.index.isin(tx.index)
    s = sf[spliced].join(tx)
    total_tx = s["NumReads"].sum()
    ri = s.loc[s["type"] == "retained_intron", "NumReads"].sum()
    n = meta["num_processed"]
    row = {
        "mode": mode,
        "run": run,
        "num_processed": n,
        "mapping_rate": meta["percent_mapped"] / 100,
        "decoy_fraction": meta["num_decoy_fragments"] / n,
        "unspliced_target_fraction": sf.loc[~spliced, "NumReads"].sum() / n,
        "transcript_fraction": total_tx / n,
        "retained_intron_share": ri / total_tx,
    }
    u = meta.get("unspliced")
    if u:
        row["usa_unspliced_fraction"] = u["unspliced_fraction"]
        row["usa_ambiguous_fraction"] = u["ambiguous_fraction"]
    return row, s


def avg_tx_length(s):
    """Per gene: TPM-weighted mean effective length and spliced read count."""
    w = s["TPM"] * s["EffectiveLength"]
    g = pd.DataFrame({"w": w, "tpm": s["TPM"], "reads": s["NumReads"], "gene": s["gene"]})
    g = g.groupby("gene").sum()
    return pd.DataFrame(
        {"avg_len": g["w"] / g["tpm"].where(g["tpm"] > 0), "reads": g["reads"]}
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out", type=Path)
    ap.add_argument("gtf", type=Path)
    ap.add_argument("--pairs", nargs="*", default=["SRR5693352,SRR5693353", "SRR5693399,SRR5693400"])
    ap.add_argument("--min-reads", type=float, default=250)
    a = ap.parse_args()

    tx = gtf_transcripts(a.gtf)
    modes = sorted(p.name for p in (a.out / "quant").iterdir())
    order = [m for m in ("gentrome", "intron", "premrna") if m in modes]
    rows, lengths = [], {}
    for mode in order:
        for q in sorted((a.out / "quant" / mode).iterdir()):
            row, s = load(a.out, mode, q.name, tx)
            rows.append(row)
            lengths[(mode, q.name)] = avg_tx_length(s)
    df = pd.DataFrame(rows)
    df.to_csv(a.out / "summary.tsv", sep="\t", index=False)

    print("### Per library\n")
    cols = ["mode", "run", "num_processed", "mapping_rate", "decoy_fraction",
            "unspliced_target_fraction", "retained_intron_share"]
    print(df[cols].to_markdown(index=False, floatfmt=".4f"))

    print("\n### Retained-intron share vs decoy share, across libraries\n")
    for mode in order:
        d = df[df["mode"] == mode]
        rho = stats.spearmanr(d["decoy_fraction"], d["retained_intron_share"])[0]
        print(f"- {mode}: Spearman rho = {rho:.2f} (n = {len(d)})")

    pairs = [p.split(",") for p in a.pairs]
    runs = [r for p in pairs for r in p]
    # One gene set for every mode: genes above --min-reads spliced reads in every
    # library under every mode, so the modes are compared on the same genes.
    common = None
    for mode in order:
        m = pd.concat({r: lengths[(mode, r)] for r in runs}, axis=1)
        ok = (m.xs("reads", axis=1, level=1) > a.min_reads).all(axis=1)
        ok &= m.xs("avg_len", axis=1, level=1).notna().all(axis=1)
        genes = set(ok[ok].index)
        common = genes if common is None else common & genes
    print("\n### Gene avgTxLength between technical replicate libraries\n")
    print(f"Genes with > {a.min_reads:g} spliced reads in every library under every mode: {len(common)}.\n")
    print("| mode | genes | SD log2 ratio, " + " | SD log2 ratio, ".join("/".join(p) for p in pairs)
          + " | r of site difference, replicate a vs b |")
    print("|---" * (3 + len(pairs)) + "|")
    for mode in order:
        m = pd.concat({r: lengths[(mode, r)] for r in runs}, axis=1)
        L = np.log2(m.xs("avg_len", axis=1, level=1).loc[sorted(common)])
        sds = [float((L[p[0]] - L[p[1]]).std()) for p in pairs]
        r = ""
        if len(pairs) == 2:
            d1 = L[pairs[0][0]] - L[pairs[1][0]]
            d2 = L[pairs[0][1]] - L[pairs[1][1]]
            r = f"{stats.pearsonr(d1, d2)[0]:.3f}"
        print(f"| {mode} | {len(L)} | " + " | ".join(f"{x:.4f}" for x in sds) + f" | {r} |")

    t = a.out / "times.tsv"
    if t.exists():
        tt = pd.read_csv(t, sep="\t")
        tt["max_rss_GiB"] = tt["max_rss_bytes"] / 2**30
        print("\n### Time and memory\n")
        print(tt[["step", "wall_s", "max_rss_GiB"]].to_markdown(index=False, floatfmt=".1f"))
        for mode in order:
            idx = a.out / f"index_{mode}"
            size = sum(f.stat().st_size for f in idx.iterdir() if f.is_file())
            print(f"- index_{mode}: {size / 2**30:.1f} GiB on disk")


if __name__ == "__main__":
    main()
