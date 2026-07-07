# BASEJ-DNAQC Containers (Open-Source Run)

This directory contains the Docker build definitions for every image used by the
**open-source** run of the `basej-dnaqc` pipeline (the default `pipeline_tool = opensource`).
The proprietary Sentieon path is intentionally excluded.

These images are **not hosted in any registry** — you build them **locally** from
the Dockerfiles here, then reference the resulting local image names in
`nextflow.config`. Each Dockerfile builds to a local tag of the form
`basejumper_<tag>` (e.g. `basejumper_ginkgo_0.3.1`).

To build every x86 image at once, see the build commands in each Dockerfile
header, or build them individually as shown in the [Building an image](#building-an-image)
section below.

## Images by pipeline stage

| Container folder | Image tag | Process | Notes |
|---|---|---|---|
| `ubuntu/` | `ubuntu_24.04_stable` | MERGE_MULTILANE_FASTQ | |
| `seqkit/` | `seqkit-2.13.0` | SEQKIT_SAMPLE | |
| `samtools/` | `samtools-1.23.1` | SAMTOOLS_SUBSAMPLE_CRAM, SAMTOOLS_MARKDUP | |
| `bwamem2_samtools/` | `bwamem2_samtools_2.2.1` | BWAMEM2_ALIGN | |
| `picard_addorreplacegroups/` | `picard_addorreplacegroups_3.1.2` | PICARD_METRICS | |
| `preseq_bam2mr/` | `preseq_bam2mr_0.1` | PRESEQ | |
| `ginkgo/` | `ginkgo_0.3.1` | GINKGO_BINUNSORT, GINKGO_SEGMENTATION_R, GINKGO_CNV_CALLER, GINKO_RDS_TO_FLAT, PARSE_RDS_CNV_METRICS | one Dockerfile bundles all Ginkgo scripts; serves all five processes |
| `ginko_parser/` | `ginko_parser_0.2.1` | GINKO_PARSE_OUTPUTS | |
| `custom_parabricks-metrics/` | `custom_parabricks-metrics_1.0.3` | GINKGO_BINS_TO_PARQUET | |
| `custom_r_qcplots/` | `custom_r_qcplots_0.2.1` | QC_PLOTS | |

## Public biocontainers (no Dockerfile needed)

These stages use unmodified, publicly available biocontainers pulled directly
from `quay.io`:

| Process | Image |
|---|---|
| FASTP_TRIM | `quay.io/biocontainers/fastp:0.20.1--h8b12597_0` |
| BAM_TO_BED | `quay.io/biocontainers/bedtools:2.28.0--hdf88d34_0` |
| MULTIQC | `quay.io/biocontainers/multiqc:1.33--pyhdfd78af_0` |

## Architecture notes

- The public release is **x86 (amd64) only**. All images build for x86; there are
  no arm variants.

## Building an image

Build each image locally from its folder. Each Dockerfile header lists its exact
build command. General form (local tags, no registry):

```bash
docker build -t basejumper_<tag> container/<folder>/
```

For example:

```bash
docker build -t basejumper_ginkgo_0.3.1 container/ginkgo/
docker build -t basejumper_bwamem2_samtools_2.2.1 container/bwamem2_samtools/
```

After building, point the `container = "..."` entries in `nextflow.config` at the
local image names you built (e.g. `container = "basejumper_ginkgo_0.3.1"`).
