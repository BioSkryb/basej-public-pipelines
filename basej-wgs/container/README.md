# BASEJ-WGS Containers (Open-Source Run)

This directory contains the Docker build definitions for every custom image used
by the **open-source** run of the `basej-wgsqc` pipeline (`--pipeline_tool opensource`).
The proprietary Sentieon path is intentionally excluded.

These images are **not hosted in any registry** — you build them **locally** from
the Dockerfiles here, then reference the resulting local image names in
`nextflow.config`. Each Dockerfile builds to a local tag of the form
`basejumper_<tag>` (e.g. `basejumper_bwamem2_samtools_picard_1.0.0`).

The public release is **x86 (amd64) only** (bwa-mem2 has no production ARM build),
so there are no arm image variants.

## Images by pipeline stage

| Container folder | Image tag | Process | Notes |
|---|---|---|---|
| `ubuntu/` | `ubuntu_24.04_stable` | MERGE_MULTILANE_FASTQ | |
| `seqkit/` | `seqkit-2.13.0` | SEQKIT_SAMPLE | |
| `samtools/` | `samtools-1.23.1` | SAMTOOLS_SUBSAMPLE_CRAM, SAMTOOLS_SUBSAMPLE_CRAM_PROPORTION | |
| `bwamem2_samtools_picard/` | `bwamem2_samtools_picard_1.0.0` | BWAMEM2_ALIGN_DEDUP_METRICS | fused align + markdup + Picard (bwa-mem2 2.2.1 + samtools 1.21 + picard 3.1.1 + python3) |
| `picard_addorreplacegroups/` | `picard_addorreplacegroups_3.1.2` | PICARD_METRICS_CRAM | open-source metrics for pre-aligned CRAM (Ultima) |
| `custom_parabricks-metrics/` | `custom_parabricks-metrics_1.0.3` | WGS_QC_METRICS_TO_PARQUET | |
| `custom_r_qcplots/` | `custom_r_qcplots_0.2.1` | WGS_QC_PLOTS | |

Six of these images (`ubuntu`, `seqkit`, `samtools`, `picard_addorreplacegroups`,
`custom_parabricks-metrics`, `custom_r_qcplots`) are identical to the images used
by the `basej-dnaqc` pipeline. Only `bwamem2_samtools_picard` is specific to
basej-wgs.

## Public biocontainers (no Dockerfile needed)

These stages use unmodified, publicly available biocontainers pulled directly
from `quay.io`:

| Process | Image |
|---|---|
| PICARD_COLLECTHSMETRICS (exome mode only) | `quay.io/biocontainers/picard:3.0.0--hdfd78af_0` |
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
docker build -t basejumper_bwamem2_samtools_picard_1.0.0 container/bwamem2_samtools_picard/
docker build -t basejumper_custom_r_qcplots_0.2.1 container/custom_r_qcplots/
```

After building, point the `container = "..."` entries in `nextflow.config` at the
local image names you built (e.g. `container = "basejumper_bwamem2_samtools_picard_1.0.0"`).
