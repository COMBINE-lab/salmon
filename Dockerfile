# salmon 2.0 (Rust) — multi-stage cargo build.
#
# This replaces the CMake/C++ Docker build (retained on the `cpp` branch for the
# C++ 1.x line). It builds the single `salmon` binary with cargo and copies it
# into a slim runtime image. All dependencies (incl. cf1-rs / piscem-rs / ksw2rs)
# resolve from crates.io, so the build needs only this repo as context.

# ---- builder ----------------------------------------------------------------
FROM rust:1.98-bookworm AS builder

WORKDIR /src
# Copy the whole workspace (Cargo.toml, crates/, .cargo/config.toml with the
# portable x86-64-v2 / NEON SIMD floor, etc.).
COPY . .

# .cargo/config.toml pins the portable SIMD floor; ksw2rs still selects
# AVX2/SSE4.1/NEON at runtime. Build only the binary crate.
RUN cargo build --release --locked --bin salmon \
    && strip target/release/salmon

# ---- runtime ----------------------------------------------------------------
FROM debian:bookworm-slim AS runtime

LABEL org.opencontainers.image.title="salmon" \
      org.opencontainers.image.description="salmon 2.0 — fast transcript quantification (Rust)" \
      org.opencontainers.image.source="https://github.com/COMBINE-lab/salmon" \
      org.opencontainers.image.licenses="BSD-3-Clause"

# The Rust binary links libc/libgcc only; the rest is for the workflow managers
# that run this image:
#   - ca-certificates: handy for users.
#   - procps: provides /bin/ps. Nextflow shells into the container and calls
#     `ps` to collect task metrics, failing the task outright with "Command 'ps'
#     required by nextflow to collect task metrics cannot be found" when it is
#     absent (#1064). It is the *only* one of Nextflow's required commands
#     missing from debian:bookworm-slim — bash, sed, grep, awk, tail, tee, date
#     and coreutils are all already present — so this single package is what
#     makes the image usable in a Nextflow pipeline.
RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates procps \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /src/target/release/salmon /usr/local/bin/salmon

# Sanity checks at build time. `ps` is checked here, not just installed above,
# so that dropping procps (or a base-image change) fails the build rather than
# silently shipping an image that breaks every Nextflow task at runtime.
RUN salmon --version \
    && command -v ps >/dev/null

ENTRYPOINT ["salmon"]
CMD ["--help"]
