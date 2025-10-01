ngs-pl
======

Synopsis
--------

**ngs-pl** [-h] [-i] [-l]

**ngs-pl** *COMMAND* [*OPTIONS*] [--debug]

Description
-----------

:program:`ngs-pl` is the front end to execute all NGS-Pipeline commands.

Options
-------

.. program:: ngs-pl

.. option::  -h, --help

   Show this help message and exit.
    
.. option:: -i

   Run in interactive mode, using IPython shell.
    
.. option:: -l, --list

   List all available commands and exit.
    
.. option:: --debug

   Enable debug mode for detailed logging.

Environment
-----------

.. envvar:: NGSENV_BASEDIR

   Path to the NGS-Pipeline base environment directory.

.. envvar:: NGS_PROMPT

   Prompt string to show in shell.

.. envvar:: NGS_PIPELINE_BASE

   Path to the NGS-Pipeline repository directory.

.. envvar:: NGS_PIPELINE_LOGLEVEL

   Set the logging level, with numeric values.

.. envvar:: NGS_PIPELINE_CMD_MODS

   Colon-separated list of additional directories to search for command
   modules.

.. envvar:: SNAKEMAKE_PROFILE

   Path to the Snakemake profile directory.

