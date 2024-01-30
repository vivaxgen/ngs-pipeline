
vivaxGEN NGS-Pipeline
=====================


The vivaxGEN NGS-Pipeline is a flexible pipeline for upstream processing (variant calling) of WGS or
targeted-sequencing data from paired-end short reads or singleton long reads NGS experiments.


Overview
--------

The vivaxGEN NGS-Pipeline was developed to cater the need to do upstream processing for
sequencing data produced within vivaxGEN project.
The pipeline uses snakemake for its workflow system and relies on micromamba to provide
its dependencies.

There are two modes of working with the pipeline:

* 2-step mode

  This mode allows the variant calling procecessing to be separated into 2 steps:

  * Mapping and indiviual sample genotyping

    This step involves preparing sample directory structures and running genotyping of each
    sample individually in parallel.

  * Joint variant calling

    This step involves running a joint-variant calling for all samples after each sample
    has been genotyped individually from the previous above step.

  The 2-step mode is suitable for processing high number of samples, especially from WGS
  experiments, where samples are processed incrementally or in batch.

* 1-step mode

  A 1-process mapping/variant-calling suitable for targeted sequencing, such as panel and
  amplicon sequencing, with smaller number of samples. 


Quick Installation
------------------

To install ngs-pipeline with all of its dependencies, run the following command on shell/terminal::

    "${SHELL}" <(curl -L https://raw.githubusercontent.com/vivaxgen/install/main/ngs-pl.sh)

When prompted for the directory to install the pipeline, enter the directory where the pipeline
will be installed.
Make sure that the installation process finish successfully.
Take a note of the directory where the pipeline is installed and the full path of its activation script.


Overview of ngs-pipeline Setting Up
-----------------------------------

Since ngs-pipeline is a variant-calling pipeline, it requires one to setup a proper configuration,
such as the reference sequence and any other settings, to work properly.


2-Step Mode Usage
---------------------



1-Step Mode Usage
---------------------


Extending ngs-pipeline
----------------------

The ngs-pipeline is developed in a way that it can be extended using Python and snakemake files.
All
