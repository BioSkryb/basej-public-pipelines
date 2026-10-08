# BaseJumper BASEJ-WGS (v2)

The BioSkryb BASEJ-WGS pipeline (manifest name `basej-wgsqc`) assesses
single-cell whole-genome or whole-exome libraries sequenced at depth. It aligns
reads, marks duplicates, and reports genome-wide (WGS) or target-panel (exome)
coverage, alignment, insert-size, GC-bias and duplication metrics. Every
biosample gets a QC score, and the run is summarized in a MultiQC report.

Version 2 replaces the BWA-MEM2 / Picard / Sentieon implementation of 1.x with
an open-source Rust toolchain. See [Migrating from 1.x](#migrating-from-1x).

# Pipeline Overview

- Concatenate multi-lane FASTQs (pipe-delimited lanes in the input CSV)
- Optional subsampling to `--max_total_reads` with **SeqKit** (off by default)
- Align, mark duplicates, coordinate-sort, index and collect QC metrics in one fused
  step: **minibwa-rs** (BWA-MEM-lineage aligner) | **dupblaster** (Picard
  MarkDuplicates algorithm) | **bam-insort** (in-memory sort + BAI) | **bskryb-qc**
  (single-pass alignment, insert-size, WGS-coverage, quality-yield and GC-bias metrics)
- Ultima CRAMs are not realigned: QC metrics come straight from the vendor CRAM and its
  duplicate flags
- Exome mode (`--mode exome`) adds **Picard CollectHsMetrics** on the panel targets and
  restricts bskryb-qc coverage to the targets
- Score every biosample (1-5 tiers plus PASS / Borderline / FAIL) and draw QC
  composition plots
- Write per-biosample Parquet tables and a **MultiQC** report

The pipeline runs on **x86_64** only. Custom images are built locally from the
Dockerfiles in [`container/`](container/). Picard and MultiQC come from public
`quay.io/biocontainers` images.

# Running Locally

## Requirements

- Java 17+ and [Nextflow](https://www.nextflow.io/) 25.10.x
- Docker (with the `buildx` plugin)
- AWS CLI, if your inputs or references are on S3

```
wget -qO- https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

## Build the container images

```
bash container/build_all_x86.sh
cat container/_build_logs/SUMMARY.log    # every line should read OK
```

This builds the `basejumper_*` images referenced in `nextflow.config` (SeqKit,
samtools, minibwa-rustqc, bskryb-rustqc, Parquet/metrics and R plotting images).
The Rust images compile from source and take a few minutes each.

## Resources

The fused alignment step requests 62 CPUs / 240 GB for WGS and 16 CPUs / 60 GB for
exome. The aligner needs about 7.2 GB for the GRCh38 index, and the in-memory sort
uses a quarter of the task memory before it spills to disk. On a smaller machine,
cap the requests:

```
--max_cpus 16 --max_memory 64.GB
```

# Reference Data

References are resolved under `--genomes_base` (default
`s3://bioskryb-shared-data`). For a local run, download the GRCh38 bundle once and
point `--genomes_base` at it. Under
`<genomes_base>/genomes/Homo_sapiens/NCBI/GRCh38/Annotation/` the pipeline reads:

| File | Used for |
|---|---|
| `GATK_bundle/Sequence/minibwa/genome.mbw`, `genome.l2b` | minibwa-rs alignment index |
| `GATK_bundle/Sequence/genome.fa`, `.fai`, `.dict` | GC bias, CRAM decoding, Picard |
| base-metrics and WGS-coverage intervals (`genomes.config`) | alignment and coverage metrics |
| `GATK_bundle/gcbias_windows/gcwin_basemetrics_w100.txt` | precomputed GC-bias windows |
| exome panel target/bait interval lists (`genomes.config`) | exome mode |

Only **GRCh38** is supported in v2.

# Input

The platform is auto-detected from the CSV columns (one platform per run).

**Illumina / Element (FASTQ)**: `biosampleName,read1,read2`

```
biosampleName,read1,read2
sample1,/data/sample1_R1.fastq.gz,/data/sample1_R2.fastq.gz
sample2,/data/s2_L001_R1.fastq.gz|/data/s2_L002_R1.fastq.gz,/data/s2_L001_R2.fastq.gz|/data/s2_L002_R2.fastq.gz
```

Separate multiple lanes with `|`. read1 and read2 must list the same number of lanes.

**Ultima (CRAM)**: `biosampleName,cram[,crai]`. If `crai` is omitted, `<cram>.crai`
is used.

# Usage

```
# Whole genome
nextflow run main.nf \
  --input_csv input.csv \
  --genomes_base /path/to/genomes \
  --outputDir results

# Exome
nextflow run main.nf \
  --input_csv input.csv \
  --genomes_base /path/to/genomes \
  --mode exome --exome_panel "xGen Exome Hyb Panel v2" \
  --outputDir results
```

| Option | Default | Description |
|---|---|---|
| `--input_csv` | (required) | Input CSV (see [Input](#input)) |
| `--outputDir` | `results` | Output directory |
| `--genomes_base` | `s3://bioskryb-shared-data` | Root of the reference bundle |
| `--genome` | `GRCh38` | Reference genome (GRCh38 only) |
| `--mode` | `wgs` | `wgs` or `exome` |
| `--exome_panel` | `xGen Exome Hyb Panel v2` | Exome panel: `xGen Exome Hyb Panel v2`, `TruSight One`, `TWIST`, `Agilent Clinical Exome`, `Twist Exome 2.0`, `xGen Pan-Cancer Hybridization Panel`, `Twist Alliance CNTG Hereditary Oncology Panel` |
| `--target_intervals` | panel default | Custom exome target interval list |
| `--skip_subsampling` | `true` | Run on all reads |
| `--max_total_reads` | `1000000000` | Subsampling cap when `--skip_subsampling false` |
| `--seqkit_sample_seed` / `--samtools_seed` | `12345` | Subsampling seeds (FASTQ / CRAM) |
| `--platform` | `ILLUMINA` | Read-group platform tag for FASTQ input |
| `--run_gcbias` | `true` | Compute GC-bias metrics |
| `--publish_bam` | `true` | Publish the deduplicated BAM + BAI |
| `--genome_territory` | GRCh38 value | Non-N reference bases for WGS coverage metrics |
| `--workspace`, `--workflow_id`, `--dataset_id` | | Labels used in output paths and tables |
| `--max_cpus`, `--max_memory` | | Per-task resource caps |

# Outputs

Under `--outputDir` (`<ws>` = `--workspace`, `<wf>` = `--workflow_id`,
`<mode>` = `wgs` or `wes`):

| Path | Contents |
|---|---|
| `workflow_outputs/<ws>/<wf>/reports/multiqc_report.html` | MultiQC report |
| `workflow_outputs/<ws>/<wf>/metrics/<mode>qc_metrics/` | All metrics and QC scores, one row per biosample |
| `workflow_outputs/<ws>/<wf>/index/per_biosample_status.csv` | QC verdict per biosample |
| `workflow_outputs/<ws>/<wf>/index/{metrics,dedup_metrics,bam}.csv` | Indexes of the per-sample files |
| `workflow_outputs/<ws>/<wf>/metrics/{rustqc,dedup}/` | Per-sample bskryb-qc and dupblaster outputs |
| `workflow_outputs/<ws>/<wf>/qc_plots/` | QC composition plot and score tables |
| `workflow_outputs/<ws>/<wf>/execution_info/tool_mqc_versions.yml` | Tool and container versions |
| `tables/<mode>qc_summary/...` | Per-biosample QC Parquet |
| `bam/<ws>/dna/tool=minibwa-rs/pipeline=wgsqc/` | Duplicate-marked BAMs + BAI (`--publish_bam`) |

# Migrating from 1.x

| | 1.x | 2.0 |
|---|---|---|
| Aligner | BWA-MEM2 (or Sentieon) | minibwa-rs |
| Duplicates | samtools markdup `-r` (removed) | dupblaster (flagged with 0x400, kept in the BAM) |
| Metrics | Picard / Sentieon | bskryb-qc single pass (Picard CollectHsMetrics in exome mode) |
| `--pipeline_tool` | `opensource` / `sentieon` | removed (single implementation) |
| `--architecture` | `arm` / `x86` | x86 only |
| Genomes | GRCh38, GRCm39 | GRCh38 |
| BAM output path | `bam/<ws>/dna/tool=bwa-mem2/pipeline=wgsqc/` | `bam/<ws>/dna/tool=minibwa-rs/pipeline=wgsqc/` |
| Reference bundle | BWA-MEM2 index | **adds** the minibwa index and GC-bias windows (see [Reference Data](#reference-data)) |

QC verdicts and the WGS/WES summary schemas matched 1.x (Sentieon) on our benchmark
biosamples. Known metric differences:

- Optical duplicates are not detected, so the optical-duplicate fields are 0.
- Insert-size metrics are empty for single-end (Ultima) data.
- Alignment-derived rates drift slightly because of the new aligner.

# Docker User / File Permissions

The custom images run as a non-root `appuser`. If tasks fail with
`touch: cannot touch '.command.trace': Permission denied`, run containers as your
own user with a `local.config`:

```groovy
docker {
    runOptions = '-u $(id -u):$(id -g)'
}
```

and add `-c local.config` to `nextflow run`.

# Testing

The nf-test cases (`wgsqc_GRCh38_test`, `wesqc_GRCh38_test`) run two biosamples each
against chr22 test references. The inputs live on the Wasabi-backed
`s3://bioskryb-public-data` bucket (contact BioSkryb support for access keys). Point
Nextflow at Wasabi in `~/.nextflow/config`:

```groovy
aws {
  region = 'us-east-1'
  client {
    endpoint = 'https://s3.us-east-1.wasabisys.com'
    s3PathStyleAccess = true
  }
}
```

Then:

```
bash container/build_all_x86.sh
nf-test test tests/main.nf.test --tag GRCh38
```

`NFTEST_GENOMES_BASE` overrides the test reference location (default
`s3://bioskryb-public-data/pipeline_resources/test_genomes`).

# Need Help?

If you need any help, please [submit a helpdesk ticket](https://bioskryb.atlassian.net/servicedesk/customer/portal/3/group/14/create/156).

# References

NOTE: Several studies have utilized BaseJumper pipelines as part of standard
quality control processes implemented through ResolveServices<sup>SM</sup>. While these
pipelines may not be explicitly cited, they are integral to the methodologies described.
