# BaseJumper BASEJ-DNA-QC (v2)

The BioSkryb BASEJ-DNA-QC pipeline evaluates the quality of single-cell DNA
libraries from low-pass sequencing (about 2M reads per sample). It reports
alignment, duplication, coverage, GC-bias, insert-size, library-complexity and
copy-number QC metrics, scores every biosample, and summarizes the run in a
MultiQC report. Use it to pick the libraries that are worth sequencing deeper.

Version 2 replaces the BWA-MEM2 / Picard / Sentieon implementation of 1.x with
an open-source Rust toolchain. See [Migrating from 1.x](#migrating-from-1x).

# Pipeline Overview

- Subsample to `--n_reads` (default 2M) with **SeqKit** (FASTQ) or **samtools** (Ultima CRAM)
- Concatenate multi-lane FASTQs (pipe-delimited lanes in the input CSV)
- Trim adapters and collect read QC with **fastp** (FASTQ)
- Align, mark duplicates and coordinate-sort in one streaming step:
  **minibwa-rs** (BWA-MEM-lineage aligner) | **dupblaster** (Picard MarkDuplicates
  algorithm) | **samtools sort**
- Ultima CRAMs are not realigned: they are decoded and duplicate-marked again with
  dupblaster after subsampling
- Collect alignment, insert-size, coverage, chrM, quality-yield and GC-bias metrics in a
  single pass with **bskryb-qc** (Rust)
- Estimate library complexity with **preseq** (`bam2mr` + `gc_extrap`)
- Call copy-number with a custom **Ginkgo** implementation (bedtools + Ginkgo R), and
  derive CNV MAPD / skew for QC scoring
- Score every biosample (1-5 tiers plus PASS / Borderline / FAIL) and draw QC
  composition, score-distribution and CNV-quadrant plots
- Write per-biosample Parquet tables and a **MultiQC** report

The pipeline runs on **x86_64** only. Custom images are built locally from the
Dockerfiles in [`container/`](container/). fastp, bedtools and MultiQC come from
public `quay.io/biocontainers` images.

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

This builds the `basejumper_*` images referenced in `nextflow.config`
(SeqKit, samtools, minibwa-rustqc, bskryb-rustqc, preseq, Ginkgo, Ginkgo parser,
Parquet/metrics and R plotting images). The Rust images compile from source and take
a few minutes each.

## Resources

The aligner loads the full GRCh38 minibwa index (about 7.2 GB) into memory, so allow
at least 16 GB of RAM. 8 CPUs and 28 GB is comfortable:

```
--max_cpus 8 --max_memory 28.GB
```

# Reference Data

References are resolved under `--genomes_base` (default
`s3://bioskryb-shared-data`). For a local run, download the GRCh38 bundle once and
point `--genomes_base` at it. The pipeline reads these from
`<genomes_base>/genomes/Homo_sapiens/NCBI/GRCh38/Annotation/GATK_bundle/`:

| File | Used for |
|---|---|
| `Sequence/minibwa/genome.mbw`, `genome.l2b` | minibwa-rs alignment index |
| `Sequence/genome.fa`, `genome.fa.fai` | GC bias and Ultima CRAM decoding |
| base-metrics intervals (`genomes.config`: `base_metrics_intervals`) | interval-restricted alignment metrics |
| `gcbias_windows/gcwin_basemetrics_w100.txt` | precomputed GC-bias windows |
| Ginkgo references for the chosen `--bin_size` / `--read_length` | CNV calling |

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

```
biosampleName,cram
sample1,/data/sample1.cram
```

# Usage

```
nextflow run main.nf \
  --input_csv input.csv \
  --genomes_base /path/to/genomes \
  --outputDir results \
  --max_cpus 8 --max_memory 28.GB
```

| Option | Default | Description |
|---|---|---|
| `--input_csv` | (required) | Input CSV (see [Input](#input)) |
| `--outputDir` | `results` | Output directory |
| `--genomes_base` | `s3://bioskryb-shared-data` | Root of the reference bundle |
| `--genome` | `GRCh38` | Reference genome (GRCh38 only) |
| `--n_reads` | `2000000` | Reads per biosample after subsampling |
| `--seqkit_sample_seed` / `--samtools_seed` | `12345` | Subsampling seeds (FASTQ / CRAM) |
| `--platform` | `ILLUMINA` | Read-group platform tag for FASTQ input |
| `--run_gcbias` | `true` | Compute GC-bias metrics |
| `--skip_cnv` | auto | Skip the Ginkgo CNV branch |
| `--bin_size` | `1000000` | Ginkgo bin size |
| `--read_length` | `50` | Read length used to select the Ginkgo reference |
| `--workspace`, `--workflow_id`, `--dataset_id` | | Labels used in output paths and tables |
| `--max_cpus`, `--max_memory` | | Per-task resource caps |

# Outputs

Under `--outputDir` (`<ws>` = `--workspace`, `<wf>` = `--workflow_id`):

| Path | Contents |
|---|---|
| `workflow_outputs/<ws>/<wf>/reports/multiqc_report.html` | MultiQC report |
| `workflow_outputs/<ws>/<wf>/metrics/qc_metrics/dnaqc_all_metrics.tsv` | All metrics, one row per biosample |
| `workflow_outputs/<ws>/<wf>/index/per_biosample_status.csv` | QC verdict per biosample |
| `workflow_outputs/<ws>/<wf>/index/{metrics,dedup_metrics,bam}.csv` | Indexes of the per-sample files |
| `workflow_outputs/<ws>/<wf>/metrics/{rustqc,dedup,preseq}/` | Per-sample bskryb-qc, dupblaster and preseq outputs |
| `workflow_outputs/<ws>/<wf>/cnv/` | Ginkgo RDS, SegCopy and CNV plots |
| `workflow_outputs/<ws>/<wf>/qc_plots/` | QC composition, score distribution and CNV-quadrant plots |
| `workflow_outputs/<ws>/<wf>/execution_info/tool_mqc_versions.yml` | Tool and container versions |
| `tables/dnaqc_summary/...` | Per-biosample QC Parquet |
| `tables/cnv_summary/...` | Per-biosample Ginkgo bin-level CNV Parquet |
| `bam/<ws>/dna/tool=minibwa-rs/pipeline=dnaqc/` | Subsampled, duplicate-marked BAMs + BAI |

# Migrating from 1.x

| | 1.x | 2.0 |
|---|---|---|
| Aligner | BWA-MEM2 (or Sentieon) | minibwa-rs |
| Duplicates | samtools markdup `-r` (removed) | dupblaster (flagged with 0x400, kept in the BAM) |
| Metrics | Picard CollectMultipleMetrics / Sentieon | bskryb-qc single pass |
| `--pipeline_tool` | `opensource` / `sentieon` | removed (single implementation) |
| `--architecture` | `arm` / `x86` | x86 only |
| Genomes | GRCh38, GRCm39, ARSUCD2 | GRCh38 |
| `groups` CSV column | used for plot grouping | not used |
| BAM output path | `bam/<ws>/dna/tool=bwa-mem2/pipeline=dnaqc/` | `bam/<ws>/dna/tool=minibwa-rs/pipeline=dnaqc/` |
| Reference bundle | BWA-MEM2 index | **adds** the minibwa index and GC-bias windows (see [Reference Data](#reference-data)) |

QC verdicts matched 1.x (Sentieon) on our benchmark biosamples. Known metric
differences:

- Optical duplicates are not detected, so the optical-duplicate fields are 0.
- Ultima `total_reads` counts primary reads only (about 12% lower than 1.x).
- Insert-size metrics are empty for single-end (Ultima) data.
- `pct_chrm` and other alignment-derived rates drift slightly because of the new
  aligner.
- `dnaqc_summary` has extra columns and a different column order. The index file
  is still named `index/metrics.csv`.

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

The nf-test case runs two chr22 samples against chr22 test references. The
inputs live on the Wasabi-backed `s3://bioskryb-public-data` bucket (contact
BioSkryb support for access keys). Point Nextflow at Wasabi in `~/.nextflow/config`:

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
`s3://bioskryb-public-data/pipeline_resources/test_genomes`). The test sets
`skip_cnv = true`, because chr22-only CNV skew is undefined.

# Need Help?

If you need any help, please [submit a helpdesk ticket](https://bioskryb.atlassian.net/servicedesk/customer/portal/3/group/14/create/156).

# References

- Chung, C., Yang, X., Hevner, R. F., et al. (2024). Cell-type-resolved mosaicism
  reveals clonal dynamics of the human forebrain. Nature, 629(8011), 384–392.
  [https://doi.org/10.1038/s41586-024-07292-5](https://doi.org/10.1038/s41586-024-07292-5)

- Zhao, Y., Luquette, L. J., Veit, A. D., et al. (2024). High-resolution detection
  of copy number alterations in single cells with HiScanner. BioRxiv.
  [https://www.biorxiv.org/content/10.1101/2024.04.26.587806v1.full](https://www.biorxiv.org/content/10.1101/2024.04.26.587806v1.full)

NOTE: Several studies have utilized BaseJumper pipelines as part of standard
quality control processes implemented through ResolveServices<sup>SM</sup>. While these
pipelines may not be explicitly cited, they are integral to the methodologies described.
