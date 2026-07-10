# BaseJumper Public Pipelines

Public [Nextflow](https://www.nextflow.io/) pipelines from BioSkryb Genomics for
single-cell genomic analysis. Each pipeline lives in its own top-level
subdirectory and is self-contained: it ships with everything needed to run
locally, using publicly available tools and container images.

The pipelines are designed to get you from input data to QC metrics and results
with minimal setup. Follow the shared prerequisites below once, then jump into
the pipeline you need.

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

## Prerequisites

These steps set up a local Ubuntu server to run any of the pipelines. You only
need to do this once per machine.

### Install Java

```
sudo apt-get install default-jdk
java -version
```

### Install AWS CLI

```
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
sudo ./aws/install
```

### Install Nextflow

```
wget -qO- https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

### Install Docker

```
# Add Docker's official GPG key:
sudo apt-get update
sudo apt-get install ca-certificates curl
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc

# Add the repository to Apt sources:
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
  $(. /etc/os-release && echo "$VERSION_CODENAME") stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update

sudo apt-get install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
```

## Container images

Each pipeline uses a mix of custom container images and public biocontainers.
The custom images are built locally from the Dockerfiles in the pipeline's
`container/` directory — see that directory's `README.md` for the full image
list and build steps. Public biocontainers (e.g. from `quay.io`) are pulled
automatically during the run. The public pipelines build for **x86 (amd64)**.

## Reference data

The pipelines read a reference-genome bundle (FASTA, indexes, intervals, and any
pipeline-specific references) from the location given by `--genomes_base`, which
defaults to the BioSkryb shared S3 path (`s3://bioskryb-shared-data`).

For a local or offline run, point `--genomes_base` at your local copy of the
reference bundle. To obtain the reference bundle and container images needed to
run these pipelines, email
[basejumper@bioskryb.com](mailto:basejumper@bioskryb.com). The directory must
contain the expected `genomes/...` layout (e.g.
`<genomes_base>/genomes/Homo_sapiens/NCBI/GRCh38/...`):

```
--genomes_base /path/to/local/genomes
```

This avoids per-run S3 access and lets a pipeline run fully offline once the
reference data and container images are in place.

## Running a pipeline

Each pipeline is run with `nextflow run main.nf` from inside its directory. A
typical run looks like:

```
cd <pipeline>
nextflow run main.nf \
  --input_csv $PWD/tests/data/inputs/<example>.csv \
  --genomes_base /path/to/local/genomes \
  --genome GRCh38 \
  --outputDir test \
  --max_cpus 8 --max_memory 30.GB
```

Common flags shared across pipelines:

- `--input_csv` — path to the input CSV (columns vary by pipeline; see the
  pipeline README)
- `--genomes_base` — reference bundle location (see [Reference data](#reference-data))
- `--genome` — reference genome (e.g. `GRCh38`)
- `--outputDir` — output directory
- `--max_cpus` / `--max_memory` — constrain resources for local runs

For pipeline-specific inputs, options, tool versions, and outputs, see the
`README.md` inside each pipeline's directory.

## Need help?

If you need any help, please email
[basejumper@bioskryb.com](mailto:basejumper@bioskryb.com).

## References

NOTE: Several studies have utilized BaseJumper pipelines as part of standard
quality control processes implemented through ResolveServices<sup>SM</sup>. While
these pipelines may not be explicitly cited, they are integral to the
methodologies described.
