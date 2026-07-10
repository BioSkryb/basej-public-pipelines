# BASEJ-LINEAGE Containers (Open-Source Run)

This directory contains the Docker build definitions for the custom images used by
the `basej-lineage` pipeline (matrix-consumer single-cell somatic lineage +
mutational-signatures). This pipeline consumes pre-built NR/NV/binary matrices and
per-sample VCFs; it uses **no proprietary tools** and needs **no license**.

These images are **not hosted in any public registry** — you build them
**locally** from the Dockerfiles here, then reference the resulting local image
names in `nextflow.config`. The convention is a local tag of the form
`basejumper_<tag>`. The public release is **x86 (amd64) only**.

To build every image at once:

```bash
bash build_all_x86.sh
```

## Images by pipeline stage

| Container folder | Image tag | Processes | Notes |
|---|---|---|---|
| `custom_snp_somatic_filter_sequoia_feb2026/` | `custom_snp_somatic_filter_sequoia_feb2026` | GENOTYPE_TABLE_FROM_ANNOTATED_VCF (via bcftools also), SEQUOIA_PHYLOGENY_{SNV,INDEL,BOTH}, SEQUOIA_VARIANT_PLACEMENT_{SNV,INDEL,BOTH}, POSTPROCESS_SEQUOIA_DRAWVAFHEAT_TREE_{SNV,INDEL,BOTH}, GATHER_LINEAGE_POSTPROCESS_ARTIFACTS, EXPORT_MUTSIG_SOMATIC_VARIANT_TABLE_FROM_BINARY_MATRIX, COMPUTE_MUTSIG_COVERAGE_FROM_VARIANT_TABLE, GATHER_MUTSIG_ARTIFACTS, PLOT_MATRIX_SIGNATURE_BARGRAPHS/*, MUSICAL_SIGNATURE_ANALYSIS, SIGDYN_SIGNATURE_ANALYSIS, MUSICATK_SIGNATURE_ANALYSIS, PLOT_COMPARATIVE_SIGNATURE_ACTIVITIES, CHECK_MATRIX_NON_EMPTY_* | Same image as `basej-somatic`; bundles the full R/awk/py script set (ggtree/treemut phylogeny, MuSiCaL, MutationalPatterns/SigDyn) in `/usr/local/bin`. |
| `sigprofiler-0.1.0/` | `sigprofiler-0.1.0` | SIGPROFILER_MATRIX_GENERATOR_FROM_MUTSIG_SOMATIC_VARIANT_TABLES, SIGPROFILER_ASSIGNMENT_* , MERGE_SIGNATURE_ACTIVITIES / MERGE_ACTIVITIES_* | python 3.10 + SigProfilerMatrixGenerator (incl. `output_directory`) + SigProfilerAssignment, with **GRCh38 baked in** (~3GB; the mutsig process passes no reference path). Recipe modelled on the committed `containers/sigprofiler_extractor` image. |

> `musicatk_enabled` is `false` by default: the feb2026 image ships R 4.2 and
> musicatk needs R >= 4.4. Keep it disabled with this image.

## Public biocontainers (no Dockerfile needed)

| Process(es) | Image |
|---|---|
| GENOTYPE_TABLE_FROM_ANNOTATED_VCF | `quay.io/biocontainers/bcftools:1.14--h88f3f91_0` |
| MULTIQC_LINEAGE | `quay.io/biocontainers/multiqc:1.33--pyhdfd78af_0` |

## Building an image

```bash
docker build -t basejumper_<tag> container/<folder>/
```

For example:

```bash
docker build -t basejumper_custom_snp_somatic_filter_sequoia_feb2026 container/custom_snp_somatic_filter_sequoia_feb2026/
docker build -t basejumper_sigprofiler-0.1.0 container/sigprofiler-0.1.0/
```

After building, point the `container = "..."` entries in `nextflow.config` at the
local image names you built (see `PUBLIC_RELEASE_INSTRUCTIONS.md`).
