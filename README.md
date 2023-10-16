# ngs-pipeline

Generic pipeline for NGS data upstream processing. For now, the pipeline can only perform variant calling and discovery, with micro-haplotype calling is in development.

Overview

This is a set of scripts and pipelines for working with NGS-based sequencing.

The scripts are mainly written in Python and the pipelines are written in Snakemake.

Processing WGS data witn ngs-pipeline generally only involves three steps:

Step 1 - preparing sample directory structures

Step 2 - running variant caller for each sample

Step 3 - running joint variant caller combining all samples


[Docs for available commands](docs/commands.rst)

