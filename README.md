# BaseJumper Public Pipelines

Public pipelines from BioSkryb Genomics for single-cell genomic analysis. Each
pipeline lives in its own top-level subdirectory and is self-contained, with its
own run instructions.

## Pipelines

| Pipeline | Directory | Description |
|---|---|---|
| basej-dnaqc | [`basej-dnaqc/`](basej-dnaqc/) | Single-cell DNA low-pass QC — subsample, align, dedup, CNV, per-biosample QC metrics and consensus scoring |
| basej-wgs | [`basej-wgs/`](basej-wgs/) | Single-cell WGS/WES QC — align, dedup, and collect whole-genome or hybrid-selection coverage metrics |
| basej-somatic | [`basej-somatic/`](basej-somatic/) | Single-cell somatic SNP/indel detection and heuristic-QC filtering with per-group variant tables |
| basej-lineage | [`basej-lineage/`](basej-lineage/) | Single-cell lineage/phylogeny reconstruction and COSMIC mutational-signature analysis from variant matrices |

## Repository layout

```
basej-public-pipelines/
├── README.md          # this file
└── <pipeline>/        # one directory per pipeline
    ├── README.md      # pipeline-specific run instructions and options
    ├── main.nf        # pipeline entrypoint
    ├── container/     # Dockerfiles for the custom images used by the pipeline
    ├── conf/          # pipeline configuration
    └── tests/         # example inputs and tests
```

## Getting started

Pick the pipeline you need from the table above and follow the `README.md` inside
its directory for setup, inputs, run commands, options, and outputs.

To obtain the reference-genome bundles and container images needed to run these
pipelines, email [basejumper@bioskryb.com](mailto:basejumper@bioskryb.com).

## Need help?

If you need any help, please email
[basejumper@bioskryb.com](mailto:basejumper@bioskryb.com).

## References

NOTE: Several studies have utilized BaseJumper pipelines as part of standard
quality control processes implemented through ResolveServices<sup>SM</sup>. While
these pipelines may not be explicitly cited, they are integral to the
methodologies described.
