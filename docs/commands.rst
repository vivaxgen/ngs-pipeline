Available Commands
==================

The ngs-pipeline command line is invoked using ``ngs-pl`` executable, i.e.::

  ngs-pl CMD [ARGS]

Fpr esample, to see the available commands, execute::

  $ ngs-pl showcmds

The available CMD are:

prepare_sample-directory
  Given a manifest file containing the sample codes and corresponding fastq files,
  the tool will prepare the directory structure for further processing.


run-sample-variant-caller
  Run the individual sample variant calling pipeline, which will go through mapping
  reads, processing mapped reads, and performing all necessary process up to getting
  GVCF files


run-joint-variant-caller
  Run the joint variant calling pipeline to process GVCF files into multi-sample VCF files.
  This command can be run after completing ``run-sample-variant-caller``


