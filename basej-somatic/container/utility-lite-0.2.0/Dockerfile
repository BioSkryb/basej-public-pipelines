# utility-lite container
# Provides: a lightweight toolbox for the bcftools/awk/python "glue" processes and
#           the python stats/parquet processes.
#           bcftools + htslib (bgzip/tabix), samtools, gawk, coreutils, and
#           python3 with pandas + pyarrow. No custom scripts are baked in — all
#           logic is inline in the process scripts of modules.nf / main.nf.
# Tag: utility-lite-0.2.0
#   (0.2.0 adds python3 + pandas + pyarrow to the 0.1.0 toolbox so this single
#    image also covers the former custom_parabricks-metrics python processes.)
# Processes: FILTER_CHOSEN_VARIANTS_BY_VEP_CHR, CONCAT_VEP_PERCHR_OUTPUTS,
#            SEQUOIA_SECOND_FILTER_MERGE, BUILD_FOCAL_PILEUP, CONCAT_FOCAL_PILEUP,
#            COLLECT_DOWNSTREAM_ARTIFACTS, EMIT_LINEAGE_INPUTS_MANIFEST,
#            COHORT_METRICS_SUMMARY,
#            VCF_TO_PARQUET_SOMATIC_VARIANTS (pandas+pyarrow),
#            SOMATIC_STATS_FROM_ANNOTATED_VCF (python stdlib)
#
# Build (x86/amd64 only):
#   docker build -t basejumper_utility-lite-0.2.0 .

FROM ubuntu:22.04
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        bcftools \
        tabix \
        samtools \
        gawk \
        coreutils \
        gzip \
        python3 \
        python3-pip \
        procps \
        ca-certificates \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# pandas + pyarrow for the VCF->parquet / somatic-stats processes (manylinux wheels).
RUN pip3 install --no-cache-dir pandas pyarrow

WORKDIR /
CMD ["/bin/bash"]
