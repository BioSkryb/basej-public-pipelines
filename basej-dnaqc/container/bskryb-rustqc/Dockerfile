# bskryb-rustqc (linux/amd64) - single-pass DNA QC metrics binary (bskryb-qc).
# Tag: bskryb-rustqc-x86_0.1.5
# Cargo.toml, Cargo.lock and src/ ship next to this Dockerfile.
# Build:
#   docker build --platform linux/amd64 -t basejumper_bskryb-rustqc-x86_0.1.5 .
FROM rust:1.88-alpine AS builder
RUN apk add --no-cache musl-dev
WORKDIR /build
COPY Cargo.toml Cargo.lock ./
COPY src ./src
# Static musl build (no glibc/htslib runtime dependencies).
RUN cargo build --release \
    && cp target/release/bskryb-qc /bskryb-qc \
    && strip /bskryb-qc

# Debian (glibc) runtime so host-mounted glibc tools (e.g. the aws CLI used for
# Nextflow S3 staging on AWS Batch) can run inside the container.
FROM debian:bookworm-slim
LABEL org.opencontainers.image.title="bskryb-rustqc" \
      org.opencontainers.image.version="0.1.5" \
      org.opencontainers.image.licenses="GPL-3.0-or-later"
# bash + procps: Nextflow runs tasks via /bin/bash and collects metrics with `ps`.
RUN apt-get update \
    && apt-get install -y --no-install-recommends bash procps ca-certificates \
    && rm -rf /var/lib/apt/lists/*
COPY --from=builder /bskryb-qc /usr/local/bin/bskryb-qc
# Run as a non-root user (Trivy DS-0002)
RUN useradd --create-home --shell /bin/bash appuser
USER appuser
# No ENTRYPOINT: Nextflow invokes /bin/bash on the task script.
