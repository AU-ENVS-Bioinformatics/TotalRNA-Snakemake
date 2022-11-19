# Snakemake workflow: TotalRNA-Snakemake

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/currocam/TotalRNA-Snakemake/workflows/Tests/badge.svg?branch=main)](https://github.com/currocam/TotalRNA-Snakemake/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake workflow for TotalRNA analysis from the Department of Environmental Science of Aarhus University. 

## Introduction

#### Overview

This pipeline manages large-scale TotalRNA meta-transcriptomic data for taxonomic analyses of SSU reads and mRNA ANALYSIS. The steps involved are:

1. Trim reads using [trim-galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/).
2. Filtering SSU and LSU reads using [sormerna](https://github.com/biocore/sortmerna) and [SILVA](https://www.arb-silva.de/).
3. Reconstructing ribosomal genes using [Metarib](https://github.com/yxxue/MetaRib).
4. Checking the quality of the ribosomal assembly using [QUAST](https://quast.sourceforge.net/).
5. Mapping RNA contigs to reads using [CoMW](https://github.com/anwarMZ/CoMW).
6. Classifying reads taxonomically using [BLAST](https://blast.ncbi.nlm.nih.gov/), [SILVA](https://www.arb-silva.de/) and [CREST](https://github.com/lanzen/CREST).
7. Assembling non-rRNA reads, filtering noncoding RNA, mapping mRNA reads to contigs and aligning contigs to [SWORD](https://academic.oup.com/bioinformatics/article/32/17/i680/2450775) using [CoMW](https://github.com/anwarMZ/CoMW).

## Getting started

#### Requirements:

It is best to pre-install Mamba before starting. All other dependencies will be installed automatically when running the pipeline for the first time.

```bash
conda activate base
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```

#### Usage
Activating conda environment:

```bash
conda activate snakemake
```

Clone this git repository to the location where you want to run your analysis:

```bash
git clone https://github.com/currocam/TotalRNA-Snakemake
cd TotalRNA-Snakemake
```

Copy raw fastq files into the ´reads´ directory. See [reads/README.md](reads/README.md) for more information. Now, we are going to rename those files and made symlinks to the `results/renamed` directory.

```bash
snakemake -n rename
snakemake -c1 rename
```

Check that all your samples are in `results/renamed`:

```bash
ls results/renamed_raw_reads/
```

Check that the pipeline will behave as expected by running a dry run and check the [configuration file](config/config.yaml) if not.

```bash
snakemake -n --use-conda
```

Finally, run the whole pipeline. A useful flag to add is `--keep-going` to prevent the pipeline to stop if an error occurs.

```bash
snakemake -c100 --use-conda --keep-going
```
