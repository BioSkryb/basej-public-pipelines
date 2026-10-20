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
| `--deepvariant_model_type` | `bioskrybv1` | Which BioSkryb DeepVariant model to run: `bioskrybv1`, `bioskrybv2` or `resolvemethylv1` |
| `--pipeline_tool` | | Pipeline tool label |

## Choosing the model

`--deepvariant_model_type` is the only switch needed to change models. It selects
the checkpoint and its population VCFs together from `conf/genomes.config`:

| Model type | Checkpoint | Tensor |
|---|---|---|
| `bioskrybv1` (default) | `bioskryb-af-20241102` (ResolveDNA) | 8-channel `[100,221,8]` |
| `bioskrybv2` | `bioskryb-rome-v2-1-1-20250625` (ResolveOME "rome") | 9-channel `[100,221,10]`, adds `identity` |
| `resolvemethylv1` | `bioskryb-methyl-af-20260413` (ResolveMethyl) | 9-channel `[100,221,9]`, same channel set as `bioskrybv2` |

All three are GRCh38 only, run under DeepVariant 1.8.0, and share the
`bioskryb-af-20241102` population VCFs (the `allele_frequency` channel).

The channel set must match the checkpoint's `example_info.json` or
`call_variants` fails on a tensor-shape mismatch. DeepVariant enforces this
itself: `make_examples` is always given `--checkpoint`, and in calling mode it
reads `example_info.json` from the checkpoint directory and derives the channels
from it, ignoring `--channel_list`. The `channel_list` recorded in
`genomes.config` is therefore documentation of what each checkpoint declares, not
a control.

Set `--deepvariant_model_type` and **not** `--deepvariant_model`. Overriding the
model path directly still calls variants — the channels follow whatever
checkpoint you point at — but it bypasses the model-to-population-VCF pairing
recorded in `genomes.config`, and it leaves the output path claiming the model
type you did not run.

Adding a future model is a `conf/genomes.config` addition (a block with `model`,
`population_vcfs` and `channel_list`) plus the `enum` in
`nextflow_schema.json`; the pipeline code does not change.

## Outputs

Outputs are written to `--outputDir`, including per-sample DeepVariant VCFs
(+ optional gVCFs), bcftools stats, and a MultiQC report.

VCFs are keyed by model, so different models can run over the same biosamples
concurrently without overwriting each other:

```
vcf/<workspace>/dna/tool=deepvariant-<version>-<model_type>/<sample>_deepvariant.vcf.gz
```

Running several in parallel is one launch per model, differing only in the model:

```bash
nextflow run main.nf --input_csv samples.csv --deepvariant_model_type bioskrybv1      -profile batch_dev
nextflow run main.nf --input_csv samples.csv --deepvariant_model_type bioskrybv2      -profile batch_dev
nextflow run main.nf --input_csv samples.csv --deepvariant_model_type resolvemethylv1 -profile batch_dev
```

Per-run artifacts (bcftools stats, MultiQC, the `index/*.csv` manifests) stay
under `workflow_outputs/<workspace>/<workflow_id>/`, which is already unique per
run because `workflow_id` comes from `TOWER_WORKFLOW_ID`. If you launch both
models outside Tower, give each run its own `--workflow_id` (or `--outputDir`)
so those do not collide.

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
