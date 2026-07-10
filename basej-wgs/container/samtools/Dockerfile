# Samtools 1.23.1 (linux/amd64). Built from source against bundled htslib.
# Mirrors the seqkit / seqkit-arm Dockerfile pattern.
#
# Build:
#   docker buildx build --platform linux/amd64 \
#     -t basejumper_samtools-1.23.1 \
#     --push container/samtools/

FROM ubuntu:24.04

ARG SAMTOOLS_VERSION=1.23.1
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates wget bzip2 xz-utils \
        gcc make autoconf perl \
        zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev \
        libncurses-dev \
    && wget -qO- "https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" \
        | tar xj -C /tmp \
    && cd /tmp/samtools-${SAMTOOLS_VERSION} \
    && ./configure --prefix=/usr/local --without-curses \
    && make -j"$(nproc)" \
    && make install \
    && cd / \
    && rm -rf /tmp/samtools-${SAMTOOLS_VERSION} \
    && apt-get purge -y \
        wget bzip2 xz-utils \
        gcc make autoconf perl \
        zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev \
        libncurses-dev \
    && apt-get install -y --no-install-recommends \
        zlib1g libbz2-1.0 liblzma5 libcurl4 libssl3t64 \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

ENTRYPOINT ["/usr/bin/env"]
CMD ["bash"]
