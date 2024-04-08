
vivaxGEN NGS-Pipeline
=====================


The vivaxGEN NGS-Pipeline is a flexible pipeline for upstream processing
(variant calling) of WGS or targeted-sequencing data from paired-end short
reads or singleton long reads NGS experiments.


Quick Installation
------------------

To install ngs-pipeline with all of its dependencies, run the following command
on shell/terminal::

    "${SHELL}" <(curl -L https://raw.githubusercontent.com/vivaxgen/install/main/ngs-pl.sh)

When prompted for the directory to install the pipeline, enter the directory
where the pipeline will be installed.
Make sure that the installation process finish successfully.
Take a note of the directory where the pipeline is installed and the full path
of its activation script.

The installation takes about 5-15 minutes (depending on internet speed), and
uses about 5-6 GB of storage.


Tutorial
--------------

A quick tutorial on setting up the environment and running the variant calling
process with the pipeline using *P vivax* data is `available here <docs/tutorial.rst>`_.


Quick Overview
--------------

The vivaxGEN NGS-Pipeline was developed to cater the need to do upstream
processing for sequencing data produced within vivaxGEN project.
The pipeline uses snakemake for its workflow system and relies on micromamba to
provide its dependencies.
The micromamba dependencies is arranged using mechanism as described in this
`repo <https://github.com/vivaxgen/install>`_.

There are two modes of working with the pipeline:

* multi-step mode

  The multi-step mode allows users to process samples incrementally from
  different sample batches in consistent ways, and is suitable for processing
  high number of samples from WGS experiments.

  The minimal steps are:

  1.  Preparing sample directory structure

      This step involves generating sample directory structures as working
      directory for the pipeline.
      Each sample will be processed in its own directory.
      The basic command for this step is::

        ngs-pl prepare-sample-directory -o OUTPUT_DIR -i MANIFEST_FILE INPUT_DIR

  2.  Running variant caller for each sample

      This step involves running variant caller for each sample in parallel.
      The result of this step would be a set of GVCF files for each sample.
      The basic command is::

        ngs-pl run-sample-variant-caller DIRECTORY_1 [DIRECTORY_2 ...]

  3.  Running joint variant caller combining all samples

      This step involves running a joint-variant calling for all samples after
      each sample has been variant-called individually from the previous above
      steps.
      The final result would be a set of VCF files containing the variants of
      all samples.
      The basic command is::

        ngs-pl run-joint-variant-caller -o OUTPUT_DIR DIRECTORY_1 [DIRECTORY_2 ...]

  There is a special command ``run-multistep-variant-caller`` that will
  execute the above commands consecutively and can be used if there is only
  one sample batch.
  However it is advised to perform the above commands manually if the sample
  set is big or there are several sample batches.

* single-step mode

  A single-step mapping/variant-calling suitable for targeted sequencing, such
  as panel and amplicon sequencing, with smaller number of samples.
  The basic command for this mode is::

    ngs-pl run-targeted-variant-caller -o OUTPUT_DIR INPUT_DIR/*.fastq.gz


See the documentation for available `commands <docs/commands.rst>`_.


Overview of ngs-pipeline Setting Up
-----------------------------------

Since ngs-pipeline is a variant-calling pipeline, it requires one to setup a
proper environment with reference sequences and any other settings before it
can be run properly.

Steps to performed in setting-up ngs-pipeline are:

1. Create a base environemnt directory

2. Generate and edit activation script in the base working directory

3. Prepare reference sequence, region files and other necessary files

4. Create a YAML-based configuration file, with proper values for each parameter

For further information about setting up the pipeline, see the `tutorial <docs/tutorial.rst>`_.


Multi-Step Mode Features
------------------------

The multi-step mode is developed to cater for incremental upstream processing
with several batches of samples which requires fully-parallelized processing
(such as WGS data) and flexible combination of configuration.

The required steps for this mode are *sample directory preparation step*
(step-1), *sample genotyping/variant-calling step* (step-2) and *joint variant-
calling step* (step-3).
Step-2 is the most resource and CPU intensive step, and probably takes almost
the majority of the processing time and storage space.

Some features of the multi-step mode are:

* Flexible configuration

  The pipeline can be configured based on different data sets, different sample
  batch, and even to individual samples. It employs cascading configuration
  feature, a mechanism where the configuration files named ``config.yaml`` are
  read, if exist, from base environment directory down to the sample directory,
  with configuration closer to the sample directory taking precedence.

* Support for incremental upstream processing
  
  Incremental data processing is very common in research fields that produces
  continuous batch of samples, such as research related to molecular
  surveillance.
  By using multiple steps in processing the data, results of any of the steps
  from previous batch can be used again with new batch of samples.

  For example, supposed there is initally a batch of samples needed to be
  processed.
  A user can run step-1, step-2 and step-3 to obtain final VCF files.
  When a new batch of samples is needed to be processed, the user will need
  to run only step-1 and step-2 to the new batch.
  Then, the user can perform step-3 by combining the results of step-2 of the
  previous batch and the new batch, to obtain the final VCF files from both
  batches.
  Hence, the user only needs to run step-2 on samples from the new batch, which
  would decrease the time and storage space needed.

* Fully-parallelized processing

  The pipeline will try to distribute the process across available cores, or
  available nodes if run under a cluster system with suitable job scheduler
  such as slurm.

* Support for troubleshooting errors

  With separate steps, any errors can be troubleshot prior to the next step,
  hence lessening the troubleshooting process.


Single-Step Mode Features
-------------------------

The single-step mode is provided for those that require simpler workflow for
upstream processing, such as panel variant calling (variant calling with
defined base positions to genotype).

The advantage of this mode is that it only requires a single command to perform
all necessary steps to obtain the final VCF files.


Extending ngs-pipeline
----------------------

The ngs-pipeline can be extended using Python and additional snakemake files.
The Python modules and the snakemake files in the ngs-pipeline can also be
imported to be used by other custom Python scripts and/or snakemake files.

To learn more about extending the pipeline or developing custom pipeline based
on ngs-pipeline, see `documentation <docs/extending.rst>`_.
