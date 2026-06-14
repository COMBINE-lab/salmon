# salmon 2.0.0 — the Rust rewrite

salmon 2.0.0 is the first release of a **from-scratch Rust reimplementation** of
salmon. It keeps the workflow you know (`salmon index` → `salmon quant` →
`quant.sf`) and the output formats your downstream tools already read, while
shipping as a single portable binary that is far easier to build, install, and
maintain. This is a new major version, so it makes a few deliberate breaking
changes; please read below and consult the documentation.

## ⚠️ Breaking change: rebuild your index

**2.0 uses a new index format and cannot read C++ (pufferfish) indices.** Rebuild
with `salmon index` from 2.0. Pointing 2.0 at a C++ index — or C++ salmon at a
2.0 index — is detected and rejected with a clear, actionable error rather than
failing cryptically.

Your **results files are unchanged**: `quant.sf`, `cmd_info.json`,
`lib_format_counts.json`, and `aux_info/meta_info.json` are drop-in, and the
inferential-replicate output (bootstrap **and** Gibbs) is written in the same
format C++ salmon used — so **tximport / tximeta / fishpond / swish keep working
without changes**.

## Highlights

- **Single portable binary.** No Boost, no CMake, no system libraries to install.
  Get it via a one-line install script, `cargo install`, conda, or Docker.
- **Same workflow and outputs**, validated against C++ salmon: selective-alignment
  quantification matches the C++ result to a per-transcript Pearson ≈ **0.999**.
- **New alignment-free `--sketch` mode** implementing an even lighter-weight variant of pseudoalignment for maximum throughput when you want it.
- **Based on sshash-rs and piscem-rs** unifying the underlying index and mapping speed, and frugal size with the piscem mapper.
- **Live progress** during mapping (a progress display on an interactive terminal; arguably nicer than the C++ progress display).
- **Now BSD-3-Clause licensed** (was GPL-3.0).
- **Feature / command line changes** Some of the least-frequently used or niche features were removed (some have stubbed command line parameters) to clean up the implementation and reduce complexity for uncommon features. Please let us know if this affects your workflow as we consider new features and what features warrant being added back.

## Installation

```sh
# prebuilt binary (Linux & macOS, x86-64 & aarch64)
curl --proto '=https' --tlsv1.2 -LsSf \
  https://github.com/COMBINE-lab/salmon/releases/latest/download/salmon-cli-installer.sh | sh

# or via cargo (Rust ≥ 1.91)
cargo install salmon-cli

# or via conda (once the 2.0 recipe lands on bioconda)
conda install -c bioconda -c conda-forge salmon

# or Docker
docker run --rm combinelab/salmon:latest salmon --version
```

Prebuilt binaries are attached to this release for `x86_64`/`aarch64` Linux and
macOS.

## Migrating from C++ salmon

See **[MIGRATION.md](https://github.com/COMBINE-lab/salmon/blob/master/MIGRATION.md)**
for the full flag-by-flag mapping. The essentials:

- **Rebuild your index** (above).
- **`salmon alevin` is removed** — single-cell quantification now lives in the
  [alevin-fry](https://github.com/COMBINE-lab/alevin-fry) ecosystem; `salmon
  alevin …` prints a redirect.
- A handful of options are **removed** (with a helpful error) or **accepted and
  ignored** where the behavior is now the default or handled differently — the
  migration guide lists each.
- `salmon swim` still swims. 🐟

## Performance & accuracy

salmon 2.0 was developed and validated by cross-checking against C++ salmon
throughout. On our benchmarks (human GEUVADIS `ERR188044` against a GRCh38 cDNA
index, and yeast `ERR458493`), selective-alignment quantification reproduces C++
salmon to per-transcript Pearson ≈ 0.999, index construction is faster, and peak
memory is lower. The new `--sketch` mode is competitive in accuracy with
dedicated pseudoaligners while integrating salmon's full bias-aware abundance
model.

## The final C++ release

salmon **1.12.0** is the **last C++ release**. Its source is preserved on the
[`cpp`](https://github.com/COMBINE-lab/salmon/tree/cpp) branch, and it remains
installable as the **`salmon-cpp`** conda package for reproducibility or
emergency fixes. New development happens here, in Rust.

## Documentation

Full docs — installation, library types, selective-alignment vs. sketch mode,
inferential replicates, the migration guide, the CLI reference, and a precise
specification of every output file — are at
**<https://combine-lab.github.io/salmon>**.

## Citation

If you use salmon, please cite:

> Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017).
> Salmon provides fast and bias-aware quantification of transcript expression.
> *Nature Methods*, 14(4), 417–419. https://doi.org/10.1038/nmeth.4197

## Acknowledgements

Thank you to the salmon user community, and to everyone whose bug reports,
feature requests, and feedback over the C++ years shaped this rewrite.
