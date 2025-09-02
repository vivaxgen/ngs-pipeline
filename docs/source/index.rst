vivaxGEN NGS-Pipeline Documentation
===================================

vivaxGEN NGS-Pipeline is an open-source, non-opinionated pipeline for variant
calling (upstream/secondary processing) of paired-end short reads or singleton
long reads NGS data.

As a non-opinionated pipeline, vivaxgen NGS-Pipeline supports the following
program/system:

* read trimmers: fastp, cutadapt, chopper

* mappers: bwa/bwa-mem2, bowtie2, minimap2

* variant callers: GATK4 (following GATK Best Practises as much as possible),
  and FreeBayes (BCFTools mpileup/call and Clair3 are currently in the works)

vivaxGEN NGS-Pipeline can be used for variant discovery and targeted/panel
variant calling.
Variant discovery is suitable for various studies such as population genetics
and GWAS, usually from WGS data.
Targeted/panel variant calling is suitable for amplicon-sequencing data
and experiments that do not require joint-variant calling such as single sample
analysis for detecting status of certain variants (reference/alternate alleles,
heterozygotes or no data).

vivaxGEN NGS-Pipeline can be installed on laptops, servers or cluster/HPC
system without the need of administrator/root privileges.
The only requirement is a UNIX-based system supported by 
`micromamba <https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html>`_
(eg. various Linux distributions, WSL2, MacOSX) with preinstalled ``curl``
and ``bash`` (these two programs are usually installed in the base system).

The pipeline can be installed in any directory.
The installation process will not clutter the home directory of the user who
performs the installation, with the exception of additional environment line in
``~/.conda/environments.txt`` (if the file already exists).
The additional line can be removed manually without affecting the working of the
pipeline.
Once installed, the pipeline can be used by any users who have read access to
the installation directory.
To uninstall the pipeline, simply remove the whole installation directory.


Quick Installation
------------------

To install ngs-pipeline with all of its dependencies, run the following command
on shell/terminal:

.. code-block:: console

    "${SHELL}" <(curl -L https://raw.githubusercontent.com/vivaxgen/ngs-pipeline/main/install.sh)

.. warning::

    If you are inside an active Conda/Mamba environemnt, please deactivate
    first before running the above command as the Conda/Mamba enviroment
    might interfere the installation process.

When prompted for the directory to install the pipeline, enter the directory
where the pipeline will be installed.
Make sure that the installation process finish successfully.
Take a note of the directory where the pipeline is installed and the full path
of its activation script.

The installation takes about 5-15 minutes (depending on internet speed), and
uses about 5-6 GB of storage.

The installation process will also detect and set the profile of the system
being installed.
If installed in HPC or cluster systems, it will try to automatically set the
pipeline to work with the correct workload manager/job scheduler (currently
only SLURM and PBSPro, however manual setting can be done for other manager).


Updating the Pipeline
---------------------

The pipeline can be updated by executing the following command in an active
environment:

.. code-block:: console

    $VVGBIN/update-box

Do note that this only updates the pipeline but does not necessarily update
the dependencies installed by ``micromamba`` and ``python pip``.
To fully update everything, a full installation needs to be performed.


.. toctree::
   :maxdepth: 2
   :caption: User Documentation

   userdocs/getting_started.rst
   userdocs/installation.rst
   userdocs/setting_configs.rst


.. toctree::
   :maxdepth: 2
   :caption: Developer Documentation

   develdocs/implementation.rst

