# minibwa-rustqc (linux/amd64) - fused aligner + dedup + sort + QC image.
# Tag: minibwa-rustqc_0.1.5
#   - minibwa-rs 0.6-r416 (MIT): BWA-MEM-lineage aligner, x86 SSE path
#   - dupblaster 0.1.1 (MIT): streaming duplicate marking
#   - bam-insort 0.1.0 (MIT): in-memory coordinate sort + BAI index
#   - bskryb-qc 0.1.5 (GPL-3.0-or-later): single-pass QC metrics
#   - samtools (Debian package)
# The bskryb-src/ and bam-insort-src/ build inputs ship next to this Dockerfile.
# Build:
#   docker build --platform linux/amd64 -t basejumper_minibwa-rustqc_0.1.5 .
FROM rust:1.90-bookworm AS builder
RUN apt-get update \
    && apt-get install -y --no-install-recommends git \
    && rm -rf /var/lib/apt/lists/*
# ---- minibwa-rs (pinned commit; bump deliberately, then re-validate) ----
WORKDIR /build/minibwa
ARG MINIBWA_RS_COMMIT=0f4b719b46c8c232ea25ea22f628e4632d037277
RUN git clone https://github.com/henriksson-lab/minibwa-rs.git . \
    && git checkout "${MINIBWA_RS_COMMIT}" \
    && sed -i 's|libsais-rs = { version = "0.2", path = "../libsais" }|libsais-rs = { version = "0.2" }|' Cargo.toml \
    && cargo build --release --features bin \
    && strip target/release/minibwa-rs
# ---- dupblaster (pinned) ----
RUN cargo install dupblaster --version 0.1.1 --root /usr/local \
    && strip /usr/local/bin/dupblaster
# ---- bskryb-qc ----
WORKDIR /build/bskryb
COPY bskryb-src/Cargo.toml bskryb-src/Cargo.lock ./
COPY bskryb-src/src ./src
RUN cargo build --release && strip target/release/bskryb-qc
# ---- bam-insort ----
WORKDIR /build/bam-insort
COPY bam-insort-src/Cargo.toml bam-insort-src/Cargo.lock ./
COPY bam-insort-src/src ./src
RUN cargo build --release && strip target/release/bam-insort

FROM debian:bookworm-slim
LABEL org.opencontainers.image.title="minibwa-rustqc" \
      org.opencontainers.image.version="0.1.5" \
      org.opencontainers.image.licenses="MIT AND GPL-3.0-or-later"
# bash + procps: Nextflow runs tasks via /bin/bash and collects metrics with `ps`.
RUN apt-get update \
    && apt-get install -y --no-install-recommends samtools ca-certificates procps bash \
    && rm -rf /var/lib/apt/lists/*
COPY --from=builder /build/minibwa/target/release/minibwa-rs /usr/local/bin/minibwa-rs
COPY --from=builder /usr/local/bin/dupblaster /usr/local/bin/dupblaster
COPY --from=builder /build/bskryb/target/release/bskryb-qc /usr/local/bin/bskryb-qc
COPY --from=builder /build/bam-insort/target/release/bam-insort /usr/local/bin/bam-insort
# Run as a non-root user (Trivy DS-0002)
RUN useradd --create-home --shell /bin/bash appuser
USER appuser
# No ENTRYPOINT: Nextflow invokes /bin/bash on the task script.
