# BASEJ-SOMATIC Containers (Open-Source Run)

This directory contains the Docker build definitions for every custom image used
by the **open-source** run of the `basej-somatic` pipeline. The open-source run
**starts from a per-sample VCF** (e.g. DNAscope output), so the proprietary
Sentieon/DNAscope image is **not** part of this directory and is not built here
(DNAscope requires a Sentieon license — see the pipeline `README.md`).

These images are **not hosted in any public registry** — you build them
**locally** from the Dockerfiles here, then reference the resulting local image
names in `nextflow.config`. Each Dockerfile header lists its exact build command;
the convention is a local tag of the form `basejumper_<tag>`
(e.g. `basejumper_custom_snp_somatic_filter_sequoia_feb2026`).

The public release is **x86 (amd64) only**; there are no arm image variants.

To build every image at once:

```bash
bash build_all_x86.sh
```

## Images by pipeline stage

| Container folder | Image tag | Processes | Notes |
|---|---|---|---|
| `custom_snp_somatic_filter_sequoia_feb2026/` | `custom_snp_somatic_filter_sequoia_feb2026` | CREATE_TAB_NVNR, SEQUOIA_BINOM_BETABINOM_TAB_NV_NR, CONCAT_FILTER_BINOM_BETABINOM_TAB_NV_NR, CUSTOM_SOMATIC_SNPINDEL_FILTERRAWTABLES, CUSTOM_CREATE_GROUP_LEVEL_TAB_DFS, SEQUOIA_SECOND_FILTER, CUSTOM_VARIANT_FILTER_PROVENANCE, CREATE_NR_NV_MATRICES_FROM_ANNOTATED_VCF, VARIANT_PREFILTER_TABLE(_PERCHR), VAF_SPLIT_VARIANTS_HEXBIN, VARIANT_FILTER_FUNNEL, COMPARE_HQSTAT_QC_REDUNDANCY, QUANTIFY_QC_FILTER_INFLUENCE, MANDATORY_VARIANTS_QC_STATUS, EXPLORE_QC_FILTER_THRESHOLDS, ANALYZE_NR_NV_PILEUP_DEPTH, PLOT_MATRIX_SCHEME_SUMMARY, CUSTOM_RSCRIPT_SOMATICSNP_FILTER_1_SAMPLELEVEL_PROCESS_PILEUP_SAMPLE_CIGAR, PLOT_GERMLINE_PREVALENCE_DISTRIBUTIONS, PLOT_ADO_GERMLINE_COMPARISON | One image bundles the full R/awk/python script set (`/usr/local/bin`); serves all sequoia/heuristic/filter/plot processes. Consolidates the former `_0.6`, `_0.7` and `feb2026` tags. |
| `variantannotation_0.4/` | `variantannotation_0.4` | PREPROCESS_VCF, BULK_GET_VARIANTS_TO_FILTER, FILTER_VEP_GERMLINE | conda env: bcftools 1.21, tabix, ensembl-vep 111, vaex/pandas |
| `utility-lite-0.2.0/` | `utility-lite-0.2.0` | FILTER_CHOSEN_VARIANTS_BY_VEP_CHR, CONCAT_VEP_PERCHR_OUTPUTS, SEQUOIA_SECOND_FILTER_MERGE, BUILD_FOCAL_PILEUP, CONCAT_FOCAL_PILEUP, COLLECT_DOWNSTREAM_ARTIFACTS, EMIT_LINEAGE_INPUTS_MANIFEST, COHORT_METRICS_SUMMARY, VCF_TO_PARQUET_SOMATIC_VARIANTS, SOMATIC_STATS_FROM_ANNOTATED_VCF | lightweight bcftools/htslib/samtools/gawk toolbox + python3 with pandas & pyarrow (0.2.0 adds pandas/pyarrow so this image also covers the former `custom_parabricks-metrics` python steps); no baked scripts |
| `r_ado_1.2/` | `r_ado_1.2` | SUMMARIZE_ADO_INTERVALS, CONCAT_SUMMARY_ADO_INTERVALS_LABELED | R + the two ADO summary scripts (consolidates the former `_1.1` / `_1.2` tags) |

## Public biocontainers (no Dockerfile needed)

These stages use unmodified, publicly available biocontainers pulled directly
from `quay.io`:

| Process(es) | Image |
|---|---|
| BCFTOOLS_STATS, CREATE_EMPTY_BULK_VARIANTS, FILTER_CHOSEN_VARIANTS_BY_BULK, GET_VARIANTS_FROM_MERGED_VCF, GET_LIST_POS_FROM_CHOSEN_VARIANTS, FILTER_DF_NV_BY_CHOSEN_VARIANTS, SUBSET/SPLIT/SORT/MERGE VEP VCF, MERGE_PROCESSED_VCF, LIST_SAMPLES_FROM_GROUP_VCF, ANNOTATE_SAMPLE_VCF, FILTER_MERGED_VCF_BY_ML_VERDICT, FILTER_VEP_GERMLINE_CHR, ANNOTATE_VCF_SCHEME_MEMBERSHIP, MERGE_ANNOTATED_SAMPLE_VCFS, SUBSET_ANNOTATED_VCF_HQSTAT_QC_DEPTH | `quay.io/biocontainers/bcftools:1.21--h8b25389_0` |
| IDENTIFY_GERMLINE_FROM_STATS, EXTRACT_GERMLINE_PREVALENCE_TABLE, SUBSET_MERGED_VCF_HIGH_CONFIDENCE_GERMLINE_FROM_STATS, CREATE_ADO_TABLE_FROM_GERMLINE_VCF, EXTRACT_NR_NV_GT_FROM_ANNOTATED_VCF | `quay.io/biocontainers/bcftools:1.14--h88f3f91_0` |
| CUSTOM_BAM_GROUP_PILEUP | `quay.io/biocontainers/samtools:1.21--h50ea8bc_0` |
| VEP_ANNOTATE | `quay.io/biocontainers/ensembl-vep:115.2--pl5321h2a3209d_1` |
| MULTIQC_SOMATIC | `quay.io/biocontainers/multiqc:1.33--pyhdfd78af_0` |

## Proprietary image (not built — VCF-start avoids it)

| Process | Image | Notes |
|---|---|---|
| SENTIEON_DNASCOPE | `…/basejumper:sentieon-202503-02` | Sentieon DNAscope germline caller. Requires a Sentieon license. The open-source run starts from a pre-computed VCF and never invokes this process, so no Dockerfile is provided and no license is needed. |

## Building an image

Build each image locally from its folder. General form (local tags, no registry):

```bash
docker build -t basejumper_<tag> container/<folder>/
```

For example:

```bash
docker build -t basejumper_custom_snp_somatic_filter_sequoia_feb2026 container/custom_snp_somatic_filter_sequoia_feb2026/
docker build -t basejumper_variantannotation_0.4 container/variantannotation_0.4/
```

After building, point the `container = "..."` entries in `nextflow.config` at the
local image names you built (see `PUBLIC_RELEASE_INSTRUCTIONS.md`).
