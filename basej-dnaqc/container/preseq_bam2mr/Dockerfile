# preseq_bam2mr container
# Provides: preseq 2.0.3 INCLUDING the `bam2mr` tool (BAM -> MR conversion) used
# by the DNA-QC PRESEQ process. NOTE: preseq 3.x renamed `bam2mr` to `to-mr`,
# which is NOT invocation-compatible (no -seg_len), so we build 2.0.3 which ships
# the real `bam2mr`. preseq 2.0.3 bundles samtools 0.1.19; gcc 7 on ubuntu:18.04
# compiles it cleanly (the -fno-common breakage only started with gcc 10).
# Also ships samtools 1.9 (built from source; ubuntu:18.04 apt only has 1.7): the
# PRESEQ step counts/filters primary reads with `samtools view` before bam2mr.
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
        liblzma-dev \
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

# samtools 1.9 (matches the samtools reported in the pipeline's software versions).
ARG SAMTOOLS_VERSION=1.9
RUN wget -q https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && \
    ./configure --prefix=/usr/local --disable-libcurl && \
    make && make install && \
    cd / && rm -rf samtools-${SAMTOOLS_VERSION} samtools-${SAMTOOLS_VERSION}.tar.bz2

WORKDIR /data

# Run as a non-root user (Trivy DS-0002)
RUN useradd --create-home --shell /bin/bash appuser && chown appuser:appuser /data
USER appuser

CMD ["/bin/bash"]
