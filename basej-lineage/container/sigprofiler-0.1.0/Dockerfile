# sigprofiler container
# Provides: python 3.10 with SigProfilerMatrixGenerator (SBS/DBS/ID matrix
#           generation, incl. output_directory) + SigProfilerAssignment (COSMIC
#           cosmic_fit), and the GRCh38 reference genome baked in (~3GB).
# Tag: sigprofiler-0.1.0
# Processes: SIGPROFILER_MATRIX_GENERATOR_FROM_MUTSIG_SOMATIC_VARIANT_TABLES,
#            SIGPROFILER_ASSIGNMENT_* , MERGE_SIGNATURE_ACTIVITIES / MERGE_ACTIVITIES_*
#
# Recipe modelled on the committed, proven containers/sigprofiler_extractor image
# (branch feature/isg/vep_module), trimmed to what the lineage pipeline uses:
# SigProfilerMatrixGenerator + SigProfilerAssignment + baked GRCh38. The mutsig
# matrix-generator process passes NO reference path and uses output_directory=,
# so a current release + in-image GRCh38 are both required.
#
# Pinned to the versions validated on 2026-07-02 (SPMG 1.3.6 / SigProfilerAssignment 1.1.3).
#
# Build (x86/amd64 only):
#   docker build -t basejumper_sigprofiler-0.1.0 .

FROM python:3.10-slim

RUN apt-get update && apt-get install -y \
        wget \
        bzip2 \
        git \
        gcc \
        g++ \
        make \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        procps \
        && apt-get clean && rm -rf /var/lib/apt/lists/*

# SigProfilerMatrixGenerator (matrix generation) + SigProfilerAssignment (cosmic_fit).
# matplotlib/seaborn/scipy back the plotting paths; pyarrow/fastparquet for TSV/parquet IO.
RUN pip install --no-cache-dir \
    SigProfilerMatrixGenerator==1.3.6 \
    SigProfilerAssignment==1.1.3 \
    matplotlib \
    seaborn \
    scipy \
    pyarrow \
    fastparquet

# Fail fast if the installed SigProfilerMatrixGenerator lacks output_directory,
# or if SigProfilerAssignment is missing (before the ~3GB genome download).
RUN python -c "import inspect; from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as m; \
assert 'output_directory' in inspect.signature(m.SigProfilerMatrixGeneratorFunc).parameters, 'SPMG lacks output_directory'; \
from SigProfilerAssignment import Analyzer; print('SigProfiler stack OK')"

# Bake the GRCh38 reference genome into the image (~3GB) — the mutsig process
# passes no reference path, so the genome must be pre-installed.
#
# On 2026-09-25 the SigProfiler FTP replaced GRCh38.tar.gz with a rebuilt genome
# (checksums only known to SPMG >= 1.4.0) and kept the previous data as
# GRCh38_Legacy.tar.gz. SPMG 1.3.6 verifies against the previous checksums, so a
# plain install('GRCh38') now fails with "GRCh38 has not been installed". Install
# the legacy archive under the GRCh38 name instead: it is the exact data this
# image was validated with. The final step fails the build unless every
# chromosome matches SPMG 1.3.6's GRCh38 checksums.
ARG SPMG_FTP=ftp://alexandrovlab-ftp.ucsd.edu/pub/tools/SigProfilerMatrixGenerator
RUN set -eux; \
    tsb_dir="$(python -c 'from SigProfilerMatrixGenerator.scripts import ref_install; print(ref_install.reference_dir().get_tsb_dir())')"; \
    mkdir -p "${tsb_dir}"; \
    wget -q -O /tmp/GRCh38_Legacy.tar.gz "${SPMG_FTP}/GRCh38_Legacy.tar.gz"; \
    tar -xzf /tmp/GRCh38_Legacy.tar.gz -C "${tsb_dir}"; \
    rm -rf "${tsb_dir}/GRCh38"; \
    mv "${tsb_dir}/GRCh38_Legacy" "${tsb_dir}/GRCh38"; \
    rm -f /tmp/GRCh38_Legacy.tar.gz; \
    python -c "from SigProfilerMatrixGenerator import install as genInstall; genInstall.install('GRCh38')"; \
    python -c "from SigProfilerMatrixGenerator.scripts.reference_genome_manager import ReferenceGenomeManager as M; m = M(); ok = m.is_genome_installed('GRCh38'); m.print_genome_checksum_verification_report('GRCh38') if not ok else None; assert ok, 'GRCh38 checksum verification failed'; print('GRCh38 (legacy data) verified')"

# Helper script(s) invoked by the lineage MERGE_ACTIVITIES_* processes (python /scripts/...).
RUN mkdir -p /scripts
COPY ./scripts /scripts
ENV PATH="/scripts:${PATH}"

WORKDIR /
CMD ["/bin/bash"]
