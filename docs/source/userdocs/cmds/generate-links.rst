ngs-pl generate-links
=====================

Synopsis
--------

**ngs-pl** [-h] [-i] [-l]

**ngs-pl generate-links** -o/--outdir *OUTDIR* [-p/--paired] [-s/--single] [-u/--underscore] [--remove-prefix] [--remove-underscore-prefix] [--use-absolute-link] *INPUTFILES*

Description
-----------

:program:`ngs-pl generate-links` generate symbolic links in batches.
This is useful to create links from FASTQ read files, in case that the link names
need some modification.

Options
-------

.. program:: ngs-pl generate-links

.. option::  -h, --help

   Show this help message and exit.
    
.. option:: -o, --outdir OUTDIR

   Output directory to create the symbolic links.
    
.. option:: --debug

   Enable debug mode for detailed logging.

