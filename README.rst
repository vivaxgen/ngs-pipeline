
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

* 3-step mode

  The 3-step mode allows users to process samples in batches or incrementally in consistent ways, which is suitable for processing high number of samples from WGS experiments.

  The steps are:

  1.  Preparing sample directory structure

      This step involves generating sample directory structures as working directory for the pipeline.
      Each sample will be processed in its own directory.
      The basic command for this step is::

        ngs-pl prepare-sample-directory -o OUTPUT_DIR -i MANIFEST_FILE INPUT_DIR

  2.  Running variant caller for each sample

      This step involves running variant caller for each sample in parallel.
      The result of this step would a set of GVCF files for each sample.
      The basic command is::

        ngs-pl run-sample-variant-caller DIRECTORY_1 [DIRECTORY_2 ...]

  3.  Running joint variant caller combining all samples

      This step involves running a joint-variant calling for all samples after each sample
      has been variant-called individually from the previous above steps.
      The final result would be a set of VCF files containing the variants of all samples.
      The basic command is::

        ngs-pl run-joint-variant-caller -o OUTPUT_DIR DIRECTORY_1 [DIRECTORY_2 ...]

* 1-step mode

  A 1-step mapping/variant-calling suitable for targeted sequencing, such as panel and
  amplicon sequencing, with smaller number of samples.
  The basic command for this mode is::

    ngs-pl run-targeted-variant-caller -o OUTPUT_DIR INPUT_DIR/*.fastq.gz

See the documentation for available `commands <docs/commands.rst>`_


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

Since ngs-pipeline is a variant-calling pipeline, it requires one to setup a proper
configuration, such as the reference sequence and any other settings, to work properly.

[TBW]


3-Step Mode Usage
---------------------

[TBW]


1-Step Mode Usage
---------------------

[TBW]


Extending ngs-pipeline
----------------------

The ngs-pipeline can be extended using Python and additional snakemake files.
The Python modules and the snakemake files in the ngs-pipeline can also be imported to be used
by other custom Python scripts and/or snakemake files.
