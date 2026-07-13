# BASEJ-GOOGLE-DEEPVARIANT Containers

This pipeline uses **only publicly available container images** — no custom
Dockerfiles need to be built.

## Images by pipeline stage

| Process | Public image | Notes |
|---|---|---|
| DEEPVARIANT_MAKE_EXAMPLES_ONLY | `google/deepvariant:1.8.0` | Official Google DeepVariant CPU image |
| DEEPVARIANT_CALL_VARIANTS | `google/deepvariant:1.8.0-gpu` | Official Google DeepVariant GPU image (requires NVIDIA GPU + CUDA) |
| DEEPVARIANT_POSTPROCESS | `google/deepvariant:1.8.0` | Same CPU image as make_examples |
| BCFTOOLS_STATS | `quay.io/biocontainers/bcftools:1.21--h8b25389_0` | Public biocontainer |
| MULTIQC_DV | `quay.io/biocontainers/multiqc:1.33--pyhdfd78af_0` | Public biocontainer |

## No build step required

All images are pulled automatically by Nextflow from their public registries:
- **Google DeepVariant**: `google/deepvariant` on Docker Hub
- **Biocontainers**: `quay.io/biocontainers/` on Quay.io

## GPU requirement

The `DEEPVARIANT_CALL_VARIANTS` process requires an NVIDIA GPU with CUDA support.
Ensure your Docker runtime is configured with GPU access (`--gpus all` is set via
`containerOptions` in `nextflow.config`).

## Architecture notes

- Google DeepVariant images are **x86 (amd64) only**.
- The GPU image requires NVIDIA drivers and the `nvidia-container-toolkit`.
