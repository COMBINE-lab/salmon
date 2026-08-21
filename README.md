# salmon

[![CI](https://github.com/COMBINE-lab/salmon/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/COMBINE-lab/salmon/actions/workflows/ci.yml)
[![GitHub release](https://img.shields.io/github/v/release/COMBINE-lab/salmon)](https://github.com/COMBINE-lab/salmon/releases/latest)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/salmon/README.html)
[![Homebrew version](https://img.shields.io/homebrew/v/salmon)](https://formulae.brew.sh/formula/salmon)
[![Docker pulls](https://img.shields.io/docker/pulls/combinelab/salmon)](https://hub.docker.com/r/combinelab/salmon)
[![License: BSD-3-Clause](https://img.shields.io/github/license/COMBINE-lab/salmon)](https://github.com/COMBINE-lab/salmon/blob/master/LICENSE)
[![DOI](https://img.shields.io/badge/DOI-10.1038%2Fnmeth.4197-blue)](https://doi.org/10.1038/nmeth.4197)

> [!IMPORTANT]
> **This is salmon 2.0 — a from-scratch Rust rewrite of salmon.** It keeps the
> same workflow (`salmon index` → `salmon quant` → `quant.sf`) and the same
> output formats downstream tools read, but it is a new major version with some
> breaking changes. The most important one: **the index format changed, so you
> must rebuild your index.** See **[MIGRATION.md](MIGRATION.md)** for the full
> list of changed/removed/new options.
>
> The final **C++** release (salmon `1.12.0`) lives on the [`cpp`](https://github.com/COMBINE-lab/salmon/tree/cpp)
> branch and remains installable as the `salmon-cpp` conda package.
>
> Single-cell quantification moved to the [alevin-fry](https://github.com/COMBINE-lab/alevin-fry)
> ecosystem (`salmon alevin` is removed).

## What is salmon?

`salmon` is a **wicked**-fast program for highly-accurate, transcript-level
quantification from RNA-seq data. It pairs a fast mapping stage — *selective
alignment*, or alignment-free *sketch* mode (`--sketch`) — with a
massively-parallel statistical model (EM/VBEM over equivalence classes) to
estimate transcript abundances. You can give salmon raw sequencing reads, or
regular alignments to the transcriptome (an **unsorted** BAM), and it uses the
same inference engine either way.

salmon 2.0 ships as a single portable binary: no compiler, Boost, or system
libraries to install.

## Installation

```sh
# install script (prebuilt binaries: Linux & macOS, x86-64 & aarch64)
curl --proto '=https' --tlsv1.2 -LsSf \
  https://github.com/COMBINE-lab/salmon/releases/latest/download/salmon-cli-installer.sh | sh

# or via cargo (Rust ≥ 1.91)
cargo install salmon-cli

# or via conda
conda install -c bioconda -c conda-forge salmon

# or Docker
docker run --rm combinelab/salmon:latest salmon --version
```

## Quick start

```sh
# 0) Unpack the sample data and go to it
tar xzf sample_data.tgz
cd sample_data

# 1) build a reusable index from a transcriptome
salmon index -t transcripts.fasta -i salmon_index -p 16

# 2) quantify (-l A auto-detects the library type)
salmon quant -i salmon_index -l A \
  -1 reads_1.fastq -2 reads_2.fastq -p 4 -o sample_quant
```

Results land in `sample_quant/quant.sf` (drop-in for tximport / tximeta /
fishpond / swish).

## Documentation

Full docs are at **<https://combine-lab.github.io/salmon>** — installation,
library types, selective-alignment vs. sketch mode, inferential replicates, the
migration guide, the CLI reference, and a precise specification of every output
file.

## Decoy-aware indexing

Accounting for [fragments of unexpected origin](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02151-8)
improves quantification. salmon can index decoy sequence (e.g. the genome)
alongside the transcriptome so reads that would otherwise be spuriously assigned
to a transcript are absorbed by the decoy. Pass a decoy-name file with
`-d/--decoys` (the decoy records must appear last in the FASTA). See the docs for
building a decoy-aware index.

## Citation

If you use salmon, please cite:

> Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017).
> Salmon provides fast and bias-aware quantification of transcript expression.
> *Nature Methods*. https://doi.org/10.1038/nmeth.4197

## License

BSD-3-Clause. See [LICENSE](LICENSE).
