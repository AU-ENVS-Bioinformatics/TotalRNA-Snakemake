# Snakemake workflow: TotalRNA-Snakemake

[![DOI](https://zenodo.org/badge/546561474.svg)](https://zenodo.org/badge/latestdoi/546561474)
[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/currocam/TotalRNA-Snakemake/workflows/Tests/badge.svg?branch=main)](https://github.com/currocam/TotalRNA-Snakemake/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake workflow for TotalRNA analysis from the Department of Environmental Science of Aarhus University. 

## TLDR

```bash
conda activate snakemake
git clone https://github.com/AU-ENVS-Bioinformatics/TotalRNA-Snakemake
cd TotalRNA-Snakemake
snakemake -c1 skip_rename # or snakemake -n rename
snakemake -c100 --use-conda --keep-going
```

## Introduction

#### Overview

This pipeline manages large-scale TotalRNA meta-transcriptomic data for taxonomic analyses of SSU reads and mRNA ANALYSIS. The steps involved are:

1. Trim reads using [trim-galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/).
2. Filtering SSU and LSU reads using [sormerna](https://github.com/biocore/sortmerna) and [SILVA](https://www.arb-silva.de/).
3. Reconstructing ribosomal genes using [Metarib](https://github.com/yxxue/MetaRib).
4. Checking the quality of the ribosomal assembly using [QUAST](https://quast.sourceforge.net/).
5. Mapping RNA contigs to reads using [BWA](https://bio-bwa.sourceforge.net/) and [samtools](https://github.com/samtools/).
6. Classifying reads taxonomically using [BLAST](https://blast.ncbi.nlm.nih.gov/), [SILVA](https://www.arb-silva.de/) and [CREST](https://github.com/lanzen/CREST).
7. Assembling non-rRNA reads ([Trinity](https://github.com/trinityrnaseq/trinityrnaseq)) and filtering noncoding RNA using the [RFam database](https://rfam.org/). 
8. Mapping mRNA contigs to reads using [BWA](https://bio-bwa.sourceforge.net/) and [samtools](https://github.com/samtools/).
9. Functional annotation of mRNA contigs using [Diamond](https://github.com/bbuchfink/diamond) and [AnnoTree](http://annotree.uwaterloo.ca/annotree/), which includes KEGG, Pfam and Tigrfam annotations for over 30,000 bacterial and 1600 archaeal genomes. 

Check the wiki for more information: https://github.com/AU-ENVS-Bioinformatics/TotalRNA-Snakemake/wiki

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

Clone this git repository to the location where you want to run your analysis. 
```bash
git clone https://github.com/AU-ENVS-Bioinformatics/TotalRNA-Snakemake TotalRNA-Snakemake-Project
cd TotalRNA-Snakemake-Project
```

Copy or symlink raw fastq files into the ´reads´ directory. See [reads/README.md](reads/README.md) for more information. Now, we are going to rename those files and made symlinks to the `results/renamed` directory. To skip this step, just copy your files into `results/renamed` and skip the next step. Alternatively, you can run `snakemake -c1 skip_rename` to symlink your files without renaming them.

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

Finally, run the whole pipeline. A useful flag to add is `--keep-going` to prevent the pipeline to stop if an error occurs. If you are running this in a shared environment, you can have all the conda environments in a shared location by adding `--conda-prefix /path/to/shared/conda/envs`. 

```bash
snakemake -c100 --use-conda --keep-going
```

## Glossary

Please find a glossary with the main generated files in the [Wiki of the project](https://github.com/AU-ENVS-Bioinformatics/TotalRNA-Snakemake).
