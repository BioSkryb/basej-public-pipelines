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
RUN python -c "from SigProfilerMatrixGenerator import install as genInstall; genInstall.install('GRCh38')"

# Helper script(s) invoked by the lineage MERGE_ACTIVITIES_* processes (python /scripts/...).
RUN mkdir -p /scripts
COPY ./scripts /scripts
ENV PATH="/scripts:${PATH}"

WORKDIR /
CMD ["/bin/bash"]
