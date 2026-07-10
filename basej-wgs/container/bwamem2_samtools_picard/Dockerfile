# BWA-MEM2 2.2.1 + samtools 1.17 + Picard 3.1.1 + python3 (linux/amd64).
#
# Fused image for the open-source basej-wgs FASTQ path
# (BWAMEM2_ALIGN_DEDUP_METRICS), which runs alignment, samtools markdup, the
# Picard metric suite and a small python3 post-processing step inside ONE
# Nextflow process so the large WGS BAM never has to be staged between tasks.
#
# bwa-mem2 has no production ARM build (x86 SSE/AVX intrinsics only), so this
# image is linux/amd64-only — the open-source path is x86-only by design.
#
# Build:
#   docker build --platform linux/amd64 -t basejumper_bwamem2_samtools_picard_1.0.0 container/bwamem2_samtools_picard/

FROM ubuntu:24.04

ARG BWAMEM2_VERSION=2.2.1
ARG DEBIAN_FRONTEND=noninteractive

# --- base tooling + BWA-MEM2 (prebuilt static binary release) ---
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates wget bzip2 \
    && wget -qO- "https://github.com/bwa-mem2/bwa-mem2/releases/download/v${BWAMEM2_VERSION}/bwa-mem2-${BWAMEM2_VERSION}_x64-linux.tar.bz2" \
        | tar xj -C /opt \
    && for b in bwa-mem2 bwa-mem2.sse41 bwa-mem2.sse42 bwa-mem2.avx bwa-mem2.avx2 bwa-mem2.avx512bw; do \
        ln -s /opt/bwa-mem2-${BWAMEM2_VERSION}_x64-linux/$b /usr/local/bin/$b; \
    done

# --- Miniforge (conda-forge default channel; no Anaconda defaults / ToS gate) ---
# Miniconda's default channels (repo.anaconda.com/pkgs/*) now require explicit
# Terms-of-Service acceptance, which breaks non-interactive builds. Miniforge
# ships conda-forge as the base channel and pulls nothing from those channels.
RUN wget --quiet "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh" -O ~/miniforge.sh \
    && /bin/bash ~/miniforge.sh -b -p /opt/conda \
    && rm ~/miniforge.sh
ENV PATH=/opt/conda/bin:$PATH

# Install only from conda-forge + bioconda, never the Anaconda defaults.
RUN conda install -y --override-channels \
        -c conda-forge -c bioconda \
        picard=3.1.1 \
        samtools=1.21 \
    && conda clean -afy

CMD ["/bin/bash"]
