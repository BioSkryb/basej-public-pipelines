# BWA-MEM2 2.2.1 + samtools 1.21 (linux/amd64).
# Bundles both tools so alignment can be piped directly into samtools sort
# within a single Nextflow process (BWAMEM2_ALIGN).
#
# Build:
#   docker buildx build --platform linux/amd64 \
#     -t basejumper_bwamem2_samtools_2.2.1 \
#     --push container/bwamem2_samtools/

FROM ubuntu:24.04

ARG BWAMEM2_VERSION=2.2.1
ARG SAMTOOLS_VERSION=1.21
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates wget bzip2 xz-utils \
        gcc g++ make autoconf perl \
        zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev \
        libncurses-dev \
    # --- BWA-MEM2 (prebuilt static binary release) ---
    && wget -qO- "https://github.com/bwa-mem2/bwa-mem2/releases/download/v${BWAMEM2_VERSION}/bwa-mem2-${BWAMEM2_VERSION}_x64-linux.tar.bz2" \
        | tar xj -C /opt \
    && ln -s /opt/bwa-mem2-${BWAMEM2_VERSION}_x64-linux/bwa-mem2 /usr/local/bin/bwa-mem2 \
    && ln -s /opt/bwa-mem2-${BWAMEM2_VERSION}_x64-linux/bwa-mem2.sse41 /usr/local/bin/bwa-mem2.sse41 \
    && ln -s /opt/bwa-mem2-${BWAMEM2_VERSION}_x64-linux/bwa-mem2.sse42 /usr/local/bin/bwa-mem2.sse42 \
    && ln -s /opt/bwa-mem2-${BWAMEM2_VERSION}_x64-linux/bwa-mem2.avx /usr/local/bin/bwa-mem2.avx \
    && ln -s /opt/bwa-mem2-${BWAMEM2_VERSION}_x64-linux/bwa-mem2.avx2 /usr/local/bin/bwa-mem2.avx2 \
    && ln -s /opt/bwa-mem2-${BWAMEM2_VERSION}_x64-linux/bwa-mem2.avx512bw /usr/local/bin/bwa-mem2.avx512bw \
    # --- samtools (built from source against bundled htslib) ---
    && wget -qO- "https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" \
        | tar xj -C /tmp \
    && cd /tmp/samtools-${SAMTOOLS_VERSION} \
    && ./configure --prefix=/usr/local \
    && make -j"$(nproc)" \
    && make install \
    && cd / \
    && rm -rf /tmp/samtools-${SAMTOOLS_VERSION} \
    # --- trim build deps, keep runtime libs ---
    && apt-get purge -y \
        gcc g++ make autoconf \
        zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev \
        libncurses-dev \
    && apt-get install -y --no-install-recommends \
        zlib1g libbz2-1.0 liblzma5 libcurl4 libssl3t64 libncurses6 \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

# Run as a non-root user (Trivy DS-0002)
RUN useradd --create-home --shell /bin/bash appuser
USER appuser

ENTRYPOINT ["/usr/bin/env"]
CMD ["bash"]
