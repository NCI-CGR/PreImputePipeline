# PreImpute Pipeline v2
# Version: 2.0.0

This repository contains a Snakemake-based preprocessing pipeline for preparing genotype data for imputation.
The pipeline performs build conversion using update-build or liftover, sample/SNP QC, HRC reference checking, and per-chromosome VCF generation.

## Features

- Supports three build modes:
  - `liftover` — lift variants from hg19 to hg38 using triple-liftOver
  - `updatebuild` — update build coordinates using a provided chain/strand file
  - `none` — pass through the input PLINK dataset unchanged
- Standard QC steps with PLINK:
  - sample call rate (`MIND`)
  - SNP call rate (`GENO`)
  - minor allele frequency (`MAF`)
  - Hardy-Weinberg equilibrium (`HWE`)
- HRC reference check using `HRC-1000G-check-bim.pl`
- Generates gzipped VCFs and `.tbi` indexes per chromosome
- Produces QC summary  manifest

## Repository structure

- `Snakefile` — main Snakemake workflow
- `config.yaml` — pipeline configuration template
- `envs/buildtools.yaml` — conda environment for build/conversion tools
- `envs/imputation.yaml` — conda environment for PLINK/BCFtools steps

## Requirements

- `snakemake`
- `conda` / `mamba`
- `perl`
- `plink` 1.9
- `bcftools`

## Required input files

- A PLINK dataset prefix with:
  - `<plink_genotype_file>.bed`
  - `<plink_genotype_file>.bim`
  - `<plink_genotype_file>.fam`
- Required scripts directory containing:
  - `triple-liftOver-main/`
  - `HRC-1000G-check-bim.pl`
  - `update_build.sh`
- Reference panel file for HRC check
- Optional chain/strand file for `updatebuild` (downloaded from https://www.strand.org.uk/)

## Configuration

Edit `config.yaml` with your dataset and pipeline settings.

Example values:

```yaml
plink_genotype_file: "/path/to/subjects"
req_scripts: "/path/to/scripts"
ref_file: "/path/to/reference_panel"
mode: "liftover"
strand_chain_file: "/path/to/file.chain"
mind: 0.05
geno: 0.05
maf: 0.01
hwe: 1e-6
include_chrX: false
include_mt: false
extra_chroms: []
pop: "ALL"
```

### Mode behavior

- `liftover`
  - performs triple-liftover from hg19 to hg38
  - requires `triple-liftOver-main` and a strand/chain file
- `updatebuild`
  - runs `update_build.sh` to update coordinates
  - requires `strand_chain_file`
- `none`
  - copies the input PLINK dataset to the pipeline build prefix without coordinate changes

## Running the pipeline

From the repository root, run:

```bash
snakemake --cores <N>
```

If you want to force a full rebuild and rerun everything:

```bash
snakemake --cores <N> --rerun-incomplete --forceall
```

## Conda environments

The pipeline uses two environment definitions:

- `envs/buildtools.yaml`
- `envs/imputation.yaml`

These are declared in the workflow using Snakemake's `conda:` directive.

Install them with Snakemake automatically or create them manually.

## Outputs

The pipeline produces:

- `run_manifest.txt` — execution summary and QC notes
- `qc_stats/*.txt` — QC reports for each filtering step
- `vcf_MIS/subjects-updated-chr<chrom>.vcf.gz` — per-chromosome gzipped VCFs
- `vcf_MIS/subjects-updated-chr<chrom>.vcf.gz.tbi` — VCF indexes

## Notes

- The pipeline automatically includes chromosomes `1..22` by default.
- Enable `include_chrX` in `config.yaml` to add chromosome X.
- Enable `include_mt` in `config.yaml` to add mitochondrial chromosome `MT`.
- Additional chromosomes can be added using `extra_chroms`.
