# BaseJumper BASEJ-LINEAGE

The BioSkryb BASEJ-LINEAGE pipeline reconstructs single-cell somatic **lineage
(phylogeny)** and **mutational signatures** from pre-built variant matrices. It is
a matrix-consumer pipeline: given paired NR/NV (variant/reference read) matrices,
a binarized variant×sample matrix, and per-sample annotated VCFs, it builds
phylogenetic trees (Sequoia/treemut), places variants on the tree, draws VAF /
digital heatmaps, and derives COSMIC mutational-signature activities
(SigProfiler, plus optional MuSiCaL and SigDyn).

This pipeline uses **only open-source tools** and public/locally-built container
images — **no license is required**. It is typically run downstream of
`basej-somatic`, which emits the NR/NV/binary matrices and annotated VCFs it
consumes.

All custom container images are **built locally** from the Dockerfiles in
[`container/`](container/README.md). Public biocontainers (bcftools, multiqc) are
pulled directly from `quay.io`.

# Pipeline Overview

Two branches run from the input matrices:

**Lineage / phylogeny**
- Build per-sample genotype tables from annotated VCFs (**bcftools**)
- Build phylogenetic trees for SNV / indel / both from NR/NV matrices (Sequoia + treemut, custom R)
- Place variants on the tree and draw VAF + digital heatmaps (custom R + poppler)
- Gather lineage figures into a report bundle

**Mutational signatures**
- Export a somatic variant table from the binary matrix
- Generate SBS/DBS/ID mutation matrices (**SigProfilerMatrixGenerator**)
- Assign COSMIC signature activities (**SigProfilerAssignment**), merge + plot
- Optional: MuSiCaL (Park lab) and SigDyn (Goncalves lab) refitting
- Aggregate into a **MultiQC** report

# Running Locally

## Install Java, AWS CLI, Nextflow, Docker

```
sudo apt-get install default-jdk
```
```
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip && sudo ./aws/install
```
```
wget -qO- https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```
Docker: see https://docs.docker.com/engine/install/ubuntu/.

## Build the container images

The custom images are built locally (x86/amd64 only) from [`container/`](container/README.md):

```
bash container/build_all_x86.sh
```

> **No license required.** This pipeline uses no proprietary tools.

# Inputs

Provide the pre-built matrices and per-sample VCFs (these are produced by the
`basej-somatic` pipeline):

- `--nr_matrix` / `--nv_matrix` — paired NR/NV matrices (rows = `CHROM_POS_REF_ALT`, cols = samples). **Required.** Use the `basej-somatic` **`PlusMandatoryNonEmpty`** scheme from the `downstream_inputs/downstream_inputs_<group>/` bundle (`<group>_NR_...PlusMandatoryNonEmpty.tsv`), **not** the `ForPhylogeny` scheme: the pre-selected `ForPhylogeny` matrix is placed in full by the phylogeny step, which leaves the variant-placement step empty and skips the VAF/digital heatmaps.
- `--binary_matrix` — binarized variant×sample matrix (0/1), same `PlusMandatoryNonEmpty` scheme. Required for the mutational-signature branch (skipped when unset).
- Per-sample annotated VCFs, via **either**:
  - `--input_csv` — CSV with header `biosampleName,vcf,vcf_index`, or
  - `--vcf_dir` — a directory of `*_somatic_annotated.vcf.gz` (+ `.tbi`) pairs.
- `--genotype_bin` (optional) — pre-binarized genotype matrix (0/0.5/1); when set, VAF discretization is skipped.

# Reference Data

Reference genome bundle is read from `--genomes_base` (default
`s3://bioskryb-shared-data`). For a local run, point it at your local copy:

```
--genomes_base /path/to/local/genomes
```

# Test Pipeline Execution

```
nextflow run main.nf \
  --nr_matrix   $PWD/tests/data/inputs/NR_matrix.tsv \
  --nv_matrix   $PWD/tests/data/inputs/NV_matrix.tsv \
  --binary_matrix $PWD/tests/data/inputs/binary_matrix.tsv \
  --input_csv   $PWD/tests/data/inputs/input_lineage.csv \
  --cohort_id   cohort \
  --gender      female \
  --genome      GRCh38 \
  --outputDir   test \
  --max_cpus 8 --max_memory 30.GB
```

> The public release is **x86 (amd64) only**.

## Command Options

```
    Usage:
        nextflow run main.nf [options]

    Script Options: see nextflow.config

        [required]
        --nr_matrix     FILE    Pre-built NR (reference read) matrix TSV
        --nv_matrix     FILE    Paired NV (variant read) matrix TSV
        --input_csv     FILE    CSV (biosampleName,vcf,vcf_index)   [or use --vcf_dir]

        [optional]
        --binary_matrix FILE    Binarized variant x sample matrix (enables mutsig branch)
        --vcf_dir       DIR     Directory of *_somatic_annotated.vcf.gz (+ .tbi)
        --genotype_bin  FILE    Pre-binarized genotype matrix (0/0.5/1)
        --cohort_id     STR     Cohort label. DEFAULT: cohort
        --gender        STR     male | female. DEFAULT: female
        --sequoia_phylogeny_mode  STR   snv,indel,both (and/or all). DEFAULT: snv,indel,both
        --musical_enabled   BOOL    Run MuSiCaL refitting. DEFAULT: true
        --sigdyn_enabled    BOOL    Run SigDyn refitting. DEFAULT: true
        --musicatk_enabled  BOOL    Run musicatk (needs R>=4.4; not in default image). DEFAULT: false
        --genome        STR     Reference genome. DEFAULT: GRCh38
        --outputDir     DIR     Path to run output directory. DEFAULT: results
        --help          BOOL    Display help message
```

## Tool versions (open-source run)

- `bcftools: 1.14`
- `MultiQC: 1.33`
- Phylogeny / signatures R+py stack: `custom_snp_somatic_filter_sequoia_feb2026`
- SigProfiler: `sigprofiler-0.1.0`

## Outputs

Outputs are written to `--outputDir`, including phylogenetic trees, variant
placements, VAF/digital heatmaps, COSMIC signature-activity tables and plots, and
a MultiQC report.

# Need Help?

If you need any help, please [submit a helpdesk ticket](https://bioskryb.atlassian.net/servicedesk/customer/portal/3/group/14/create/156).

# References

NOTE: Several studies have utilized BaseJumper pipelines as part of standard
quality control processes implemented through ResolveServices<sup>SM</sup>. While
these pipelines may not be explicitly cited, they are integral to the
methodologies described.
