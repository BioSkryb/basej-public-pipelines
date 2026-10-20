# BaseJumper Public Pipelines

Public pipelines from BioSkryb Genomics for single-cell genomic analysis. Each
pipeline lives in its own top-level subdirectory and is self-contained, with its
own run instructions.

## Pipelines

| Pipeline | Directory | Description |
|---|---|---|
| basej-dnaqc | [`basej-dnaqc/`](basej-dnaqc/) | Single-cell DNA low-pass QC (v2, open-source Rust toolchain) — subsample, align, dedup, CNV, per-biosample QC metrics and consensus scoring |
| basej-rnaqc | [`basej-rnaqc/`](basej-rnaqc/) | Single-cell RNA library QC — STAR 2-pass alignment, HTSeq gene quantification, Qualimap genomic composition, and composite QC scoring |
| basej-wgs | [`basej-wgs/`](basej-wgs/) | Single-cell WGS/WES QC (v2, open-source Rust toolchain) — align, dedup, and collect whole-genome or hybrid-selection coverage metrics |
| basej-google-deepvariant | [`basej-google-deepvariant/`](basej-google-deepvariant/) | Germline SNV/indel calling from BAM/CRAM with Google DeepVariant using a BioSkryb custom-trained model that corrects PTA artifacts |
| basej-somatic | [`basej-somatic/`](basej-somatic/) | Single-cell somatic SNP/indel detection and heuristic-QC filtering with per-group variant tables |
| basej-lineage | [`basej-lineage/`](basej-lineage/) | Single-cell lineage/phylogeny reconstruction and COSMIC mutational-signature analysis from variant matrices |

Version 2 of basej-dnaqc and basej-wgs replaces the BWA-MEM2 / Picard / Sentieon
implementation with minibwa-rs, dupblaster and bskryb-qc, and is x86_64 / GRCh38 only.
Each pipeline README has a "Migrating from 1.x" section.

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

To obtain the reference-genome bundles and pipeline resources needed to run these
pipelines, email [basejumper@bioskryb.com](mailto:basejumper@bioskryb.com).

## Need help?

If you need any help, please email
[basejumper@bioskryb.com](mailto:basejumper@bioskryb.com).

## References

NOTE: Several studies have utilized BaseJumper pipelines as part of standard
quality control processes implemented through ResolveServices<sup>SM</sup>. While
these pipelines may not be explicitly cited, they are integral to the
methodologies described.
