# preseq_bam2mr container
# Provides: preseq 2.0.3 INCLUDING the `bam2mr` tool (BAM -> MR conversion) used
# by the DNA-QC PRESEQ process. NOTE: preseq 3.x renamed `bam2mr` to `to-mr`,
# which is NOT invocation-compatible (no -seg_len), so we build 2.0.3 which ships
# the real `bam2mr`. preseq 2.0.3 bundles samtools 0.1.19; gcc 7 on ubuntu:18.04
# compiles it cleanly (the -fno-common breakage only started with gcc 10).
# Tag: preseq_bam2mr_0.1
#
# Build:
#   docker build -t basejumper_preseq_bam2mr_0.1 .

FROM ubuntu:18.04

ENV DEBIAN_FRONTEND=noninteractive

# Build dependencies: GSL (preseq), zlib + ncurses (bundled samtools 0.1.19)
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        wget \
        ca-certificates \
        zlib1g-dev \
        libgsl-dev \
        libncurses5-dev \
        libbz2-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Build preseq 2.0.3 + bam2mr. `make all` builds both `preseq` and `bam2mr`
# (SAMTOOLS_DIR defaults to the bundled samtools/ in the tarball).
ARG PRESEQ_VERSION=2.0.3
RUN wget -q https://github.com/smithlabcode/preseq/releases/download/v${PRESEQ_VERSION}/preseq_v${PRESEQ_VERSION}.tar.bz2 && \
    tar xjf preseq_v${PRESEQ_VERSION}.tar.bz2 && \
    cd preseq && \
    make all && \
    install -m 755 preseq bam2mr /usr/local/bin/ && \
    cd / && rm -rf preseq preseq_v${PRESEQ_VERSION}.tar.bz2

WORKDIR /data

# Run as a non-root user (Trivy DS-0002)
RUN useradd --create-home --shell /bin/bash appuser && chown appuser:appuser /data
USER appuser

CMD ["/bin/bash"]
