Settings and Configuration
==========================


Environment Variables
---------------------

The activation script of NGS-Pipeline sets several environement variables
in the activated shell as following:

VVG_BASEDIR
VVGBIN
NGS_PIPELINE_BASE
NGSENV_BASEDIR

Users can set some environment variables either in a profile file under
``$VVG_BASEDIR/etc/bashrc.d/`` directory, in ``$NGSENV_BASEDIR/profile``, or
directly from the shell using notation ``ENVVAR_NAME=value``.

The following are the list of the environment variables:

.. list-table:: Environment Variables
    :header-rows: 1

    * - Variable name
      - Remarks
    * - NGS_PIPELINE_LOGLEVEL
      -
    * - NGS_PIPELINE_LOGFILE
      -
    * - NGS_PIPELINE_JOBS
      -
    * - NGS_PIPELINE_FORCE
      -
    * - NGS_PIPELINE_NO_CONFIG_CASCADE
      -
    * - NGS_PIPELINE_CMD_MODS
      -
    * - SNAKEMAKE_CLUSTER_EXTRA_FLAGS
      -


Configuration Files
-------------------

User can set custom configuration in YAML-format config files, usually named
as ``config.yaml``.


Cascading Configuration
~~~~~~~~~~~~~~~~~~~~~~~

In default mode, NGS-Pipeline will perform cascading configuration setting.
Essentially, when a snakemake workflow process is being executed, the process
will read for a ``config.yaml`` file in current working directory if the file
exists.
The workflow process will try to read another ``config.yaml`` file in the
parent directory, and the parent of the parent directory, and so on until it
reaches the NGSENV_BASEDIR as the root of the hieararchy.
Settings from ``config.yaml`` in the directory farther to the current working
directory will be overridden by the settings from ``config.yaml`` in the
directory closer to the working directory (up to the working directory itself).
With this scheme, it is easy to arrange configurations applied to whole
project, but then customized for certain sample sets down to individual sample.

To illustrate how the cascading configuration works, assume that we have the
following directory layout:

.. code-block:: console

    NGSENV_BASEDIR
    ├── config.yaml (1)
    └── sets
        ├── clinical-samples
        │   ├── analysis
        │   │   └── samples
        │   │       ├── patient-001
        │   │       │   └── config.yaml (2)
        │   │       ├── patient-002
        │   │       └── patient-003
        │   └── config.yaml (3)
        ├── joint-varcall
        │   └── config.yaml (4)
        ├── public-samples
        │   └── analysis
        │       └── samples
        │           ├── P0001
        │           └── P0002
        └── study-A
            ├── analysis
            │   └── samples
            │       ├── A001
            │       └── A002
            └── config.yaml (5)

When ``run-sample-variant-caller`` command is executed to perform per sample
processing (mapping, genotyping) to samples in the ``clinical-samples``,
``public-samples`` and ``study-A`` directories, a snakemake workflow is being
run for each sample with the respective sample directory as working directory.
For ``patient-001`` sample, the workflow will encounter ``config.yaml (2)``,
then ``config.yaml (3)``, and the base ``config.yaml (1)``.
The settings in the ``config.yaml (1)`` will be overridden with any settings in
``config.yaml (3)``, which then will be overridden by any settings in
``config.yaml (2)``.

Likewise, for sample ``A001``, the applied settings will be those from
``config.yaml (1)`` which will then be overriden by ``config.yaml (5)``.
For sample ``P0001``, the settings will only use the ones from
``config.yaml (1)``.

With this scheme, it is easy to setup general settings for all samples in
``config.yaml (1)``, set some custom settings for in-house-sequenced
``study-A`` samples (such as keeping the proper-paired bam files for SRA
submission) in ``config.yaml (5)``, set some custom settings for all clinical
samples (such as keeping the final bam files for manual inspection) in
``config.yaml (3)`` and set specific settings for just sample ``patient-001``
(such as lowering some thresholds as the sample is of lower quality) in
``config.yaml (2)``.

The cascading configuration can be opted out by using ``--no-config-cascade``
argument in most of NGS-Pipeline commands.


Configurations
~~~~~~~~~~~~~~

.. list-table:: Configurations to select workflows
    :header-rows: 1

    * - Config Name
      - Remarks
      - Default Value
      - Available Values
    * - read_trimmer_wf
      -
      - ssf_trimmer_null.smk
      - ssf_trimmer_fastp.smk ssf_trimmer_cutadapt.smk
    * - reads_mapper_wf
      -
      - ssf_mapper_bwa.smk
      - ssf_mapper_minimap2.smk ssf_mapper_bowtie2.smk
    * - variant_caller_wf
      -
      - ssf_varcall_gatk.smk
      -
    * - joint_variant_caller_wf
      -
      - jointvarcall_gatk.smk
      - jointvarcall_freebayes.smk jointvarcall_bcftools.smk jointvarcall_clair3.smk

.. list-table:: Configuration for map processing
    :header-rows: 1

    * - Config Name
      - Remarks
      - Default Value
      - Available Values
    * - refseq_file
      -
      -
      -
    * - refmap_file
      -
      -
      -
    * - deduplicate
      -
      - True
      - False
    * - keep_paired_bam
      -
      - False
      - True
    * - keep_final_bam
      -
      - False
      - True
