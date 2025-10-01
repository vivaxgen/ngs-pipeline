ngs-pl run-targeted-variant-caller
==================================

Synopsis
--------

**ngs-pl run-targeted-variant-caller** -o/--outdir *OUTDIR* [--snakefile] [--panel] [-t/--target] [-p/--paired] [-s/--single] [-u/--underscore] *INPUTFILES*

Description
-----------

:program:`ngs-pl run-targeted-variant-caller` runs targeted variant calling workdlow.


Options
-------

.. program:: ngs-pl run-targeted-variant-caller

.. option::  -h, --help

   Show this help message and exit.
    
.. option:: -o, --outdir OUTDIR

   Output directory to create the symbolic links.
    
.. option:: --debug

   Enable debug mode for detailed logging.


Configurations
--------------

.. option:: targeted_variant_caller_wf

   Snakemake workflow file (snakefile) to run.
