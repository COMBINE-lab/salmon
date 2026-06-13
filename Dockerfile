# salmon 2.0 (Rust) — multi-stage cargo build.
#
# This replaces the CMake/C++ Docker build (docker/Dockerfile, retained on the
# `cpp` branch for the C++ 1.x line). It builds the single `salmon` binary with
# cargo and copies it into a slim runtime image.
#
# NOTE: this builds once the salmon workspace depends on the *published*
# crates.io versions of cf1-rs / piscem-rs / ksw2rs (release plan R3). Before
# that, those are local path deps (../cf1-rs, ../piscem-rs, ../../../ksw2rs) and
# a build from this context alone cannot resolve them.

# ---- builder ----------------------------------------------------------------
FROM rust:1.91-bookworm AS builder

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

# The Rust binary links libc/libgcc only; ca-certificates is handy for users.
RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /src/target/release/salmon /usr/local/bin/salmon

# Sanity check at build time.
RUN salmon --version

ENTRYPOINT ["salmon"]
CMD ["--help"]
