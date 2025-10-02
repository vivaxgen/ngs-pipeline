ngs-pl run-targeted-variant-caller
==================================

Synopsis
--------

**ngs-pl run-targeted-variant-caller** -o/--outdir *OUTDIR* [--snakefile] [--panel] [-t/--target] [-p/--paired] [-s/--single] [-u/--underscore] *INPUTFILES*

Description
-----------

:program:`ngs-pl run-targeted-variant-caller` runs targeted variant calling
workflow to produce single-sample VCF file for each sample.
Each VCF file will then be merged into a multi-sample VCF file when ``--target``
is set to ``merged_vcf``.


Options
-------

.. program:: ngs-pl run-targeted-variant-caller

.. option::  -h, --help

   Show this help message and exit.
    
.. option:: -o, --outdir OUTDIR

   Output directory to create the symbolic links.
    
.. option:: --debug

   Enable debug mode for detailed logging.

.. option:: --snakefile

   Snakefile to use for the workflow.
   When supplied without absolute or relative path, the snakefile will be
   looked up from the available snakefiles included with the pipeline.
   See the configuration section below for the list of available snakefiles.

.. option:: -t, --target

   Set the target of the workflow.
   Available targets depend on the used snakefile, see the configuration
   section below.
   

Configurations
--------------

.. option:: targeted_variant_caller_wf

   Snakemake workflow file (snakefile) to run.

   Available snakefiles are:

   msf_panel_varcall_lr.smk
   
     This is the workflow for panel variant calling using singleton long reads
     from ONT or PacBio.

     Available targets are:

     merged_vcf
       This results in a multi-sample VCF file as VCF files from each sample is
       merged using ``bcftools merge`` command.

     merged_report
       This results in a Excel file containing report of genotypes of each
       sample.

   msf_panel_varcall_pe.smk

     This is the workflow panel variant calling using paired-end short reads from
     Illumina platforms.

     Available targets are the same as msf_panel_varcall_lr.smk above

.. option:: target

   The targets will be specific to each of snakefile.
   See the snakefile above to get the available targets.


