# BaseJumper Public Pipelines

Public-ready [Nextflow](https://www.nextflow.io/) pipelines from BioSkryb Genomics.
Each pipeline lives in its own top-level subdirectory and is a self-contained,
open-source (no Sentieon license required) transform of its internal source.

> **Status: PRIVATE.** This repository is temporarily private while pipelines are
> vetted and reviewed. It will be made public manually, via a console change,
> once review is complete.

## Layout

```
basej-public-pipelines/
├── README.md          # this file
└── <pipeline>/        # one directory per published pipeline
```

| Pipeline | Directory | Description |
|---|---|---|
| basej-dnaqc | [`basej-dnaqc/`](basej-dnaqc/) | Single-cell DNA QC (subsample, align, dedup, CNV, consensus) |

## How pipelines get here

Each pipeline is published by an automated, approval-gated GitHub Actions workflow
in the private `nf-bioskryb-utils` monorepo. On release, the workflow regenerates
the pipeline's public-ready form (local container tags, no internal AWS Batch /
resource-label directives, self-contained `conf/`), runs security gates
(secret scan, Dockerfile lint, dependency + image CVE scan, SBOMs), and — after a
manual approval — replaces **only that pipeline's subdirectory** here. Sibling
pipeline directories are never touched.

## Running a pipeline

See the `README.md` inside each pipeline's directory for run instructions and the
list of container images to build locally.
