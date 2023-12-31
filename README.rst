vivaxGEN NGS-Pipeline
=====================


The viaxGEN NGS-Pipeline is a flexible pipeline for upstream processing (variant calling) of WGS or
targeted-sequencing data from paired-end short reads or singleton long reads NGS experiments.


Overview
--------

The vivaxGEN NGS-Pipeline was developed to cater the need to do upstream processing for
sequencing data produced within vivaxGEN project.
The pipeline uses snakemake for its workflow system and relies on micromamba to provide
its dependencies.
There are two modes of working with the pipeline:

* a 2-process mechanism with different steps for mapping and joint variant calling which is
  suitable for whole-genome sequencing or variant discovery with big number of samples, and

* a 1-process mapping/variant-calling suitable for targeted sequencing, such
  as panel and amplicon sequencing, with smaller number of samples. 


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


2-Step WGS Mode Usage
---------------------



1-Step Targeted Mode Usage
--------------------------


Extending ngs-pipeline
----------------------

The ngs-pipeline is developed in a way that it can be extended using Python and snakemake files.
All
