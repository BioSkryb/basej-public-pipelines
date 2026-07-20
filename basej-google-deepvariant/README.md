# BaseJumper BASEJ-GOOGLE-DEEPVARIANT

The BioSkryb BASEJ-GOOGLE-DEEPVARIANT pipeline performs germline SNV/indel
calling from aligned BAM/CRAM files using Google DeepVariant with a BioSkryb
custom-trained model that corrects PTA (Primary Template-directed Amplification)
artifacts.

## Pipeline Overview

The pipeline runs the three DeepVariant stages:

1. **make_examples** — Converts aligned reads into pileup image tensors (tfrecords)
2. **call_variants** — Runs the trained model on the tfrecords to produce variant calls
3. **postprocess_variants** — Assembles the final VCF (+ optional gVCF) from model output

Additionally:
- **bcftools stats** — Collects variant-level summary statistics
- **MultiQC** — Generates a summary report

## Input

A CSV file with columns:
- `biosampleName` — sample identifier
- `bam` or `cram` — path to the aligned BAM or CRAM file

```csv
biosampleName,bam
sample1,s3://bucket/path/to/sample1.bam
sample2,s3://bucket/path/to/sample2.bam
```

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input_csv` | (required) | Input CSV with BAM/CRAM paths |
| `--genome` | `GRCh38` | Reference genome |
| `--mode` | `wgs` | `wgs` or `exome` (exome requires `--regions`) |
| `--regions` | | BED file for exome calling intervals |
| `--outputDir` | `results` | Output directory |
| `--pipeline_tool` | | Pipeline tool label |

## Outputs

Outputs are written to `--outputDir`, including per-sample DeepVariant VCFs
(+ optional gVCFs), bcftools stats, and a MultiQC report.

# Testing

## Test Data Access

Test data is stored on Wasabi-backed S3 at `s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/`.

To access the test data:

**Step 1 — Get your access keys**

Retrieve your AWS credentials from BioSkryb support (contact basejumper-support for the access link).

**Step 2 — Set environment variables**

```bash
export AWS_ACCESS_KEY_ID=<provided_access_key>
export AWS_SECRET_ACCESS_KEY=<provided_secret_key>
export AWS_DEFAULT_REGION=us-east-1
```

## Running a Test

Run the pipeline with the provided test input CSV (1 sample, ~1M reads WGS BAM, ~1 hour):

```bash
nextflow run main.nf \
  --input_csv tests/data/inputs/nftest_input.csv \
  --max_cpus 8 --max_memory 24.GB --architecture x86 \
  --genome GRCh38 --mode wgs \
  --outputDir results_test
```

## nf-test (Automated Testing)

Install nf-test (requires Java 11+):

```bash
curl -fsSL https://code.askimed.com/install/nf-test | bash
mv nf-test /usr/local/bin/
```

Run the automated tests:

```bash
# Run all tests
nf-test test

# Run only the GRCh38 test
nf-test test tests/main.nf.test --tag GRCh38
```


# Need Help?

If you need any help, please [submit a helpdesk ticket](https://bioskryb.atlassian.net/servicedesk/customer/portal/3/group/14/create/156).
