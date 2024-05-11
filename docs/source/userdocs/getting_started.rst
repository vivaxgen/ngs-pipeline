Getting Started
===============

The document will work through setting up a pipeline to process *Plasmodium
vivax* sequence data using PvP01_v1 reference sequence and perform some
processing of *P vivax* sequence data.

.. note::

  Run this tutorial under a terminal multiplexer such as ``tmux`` or ``screen``
  when using a remote system to avoid termination of a running program if the
  network is disrupted.


Quick Installation
------------------

If vivaxGEN NGS-Pipeline has not been installed in your system, run the
following command to install it:

.. prompt:: bash

    "${SHELL}" <(curl -L https://raw.githubusercontent.com/vivaxgen/ngs-pipeline/main/install.sh)

Enter the target directory at the prompt (or just press Enter to use the
default).
The installation will take about 5-15 minutes depending on the internet speed
and will use about 5-6 GB of storage.
Take a note of the directory where the pipeline is installed and the full path
of its activation script.

Try to activate the NGS-Pipeline environment by executing its activation script:

.. code-block:: console

  YOUR_INSTALATION_DIRECTORY/bin/activate

.. note::

  When an environment is activated, the shell prompt will change to reflect
  the environment name.


Checking and Setting the Profile
--------------------------------

The installation process checks and sets the correct cluster profile for
Snakemake workflow if a workload manager/job scheduler is installed in
a HPC/cluster system.

If SLURM or PBSPro is installed in the system, the installation script will
try to set up a default setting.
Some HPC settings require additional flags or arguments for job submission, in
which case the flags/arguments can be supplied with environment variable
``SNAKEMAKE_CLUSTER_EXTRA_FLAGS``, by following instructions in
:ref:`cluster_extra_flags`.
If the installer can not detect existing workload manager, please follow
:ref:`profile_manual_setup` to manually set the cluster profile.

If the installer could not detect any workload managers/job schedulers, it will
set the profiler based on the available cores and memory of the system.

To see which profile set by the installer, activate the NGS-Pipeline environment
and check the link to the Snakefile profile by running the following command:

.. code-block:: console

  ls -l $VVG_BASEDIR/etc/bashrc.d/99-snakemake-profile

The link should point to the correct profile.
If the link is not correct, adjust the link accordingly following the
instruction in :ref:`snakemake_profile_setting`.


Preparing Base Environment Directory
------------------------------------

For normal workflows, vivaxGEN NGS-Pipeline requires that all necessary files,
with the exception of the pipeline and all supporting software themselves, to
reside in a single hierarchical directory (unless the pipeline has been set up
with specific settings).
The base environment directory can be prepared with automatic method (with
preset settings) or manual method (with custom settings).

Regardless of whether using automatic or manual method, the following are
the mandatory steps to prepare the base environment directory:

#.  Activate the vivaxGEN NGS-Pipeline environment by running its activation
    script, as noted after the automatic installation finished, eg:

    .. prompt:: bash

          NGS-PIPELINE_INSTALL_DIR/activate

#.  Setup the base working directory, eg: ``/data/Pv-wgs/PvP01_v1``:

    .. prompt:: bash

          ngs-pl setup-base-directory /data/Pv-wgs/PvP01_v1
          cd /data/Pv-wgs/PvP01_v1

    .. tip::

      To easily identify and differentiate between several base environemnt
      directory, it is recommended to use the reference genome name as part
      of the directory name.

#.  Exit the current environment and activate the new environment using
    the new ``activate`` file:

    .. prompt:: bash

          exit
          /data/Pv-wgs/PvP01_v1/activate

    Once activated, the environment directory can be accessed using environment
    variable ``NGSENV_BASEDIR``.

To continue preparing the base enviroment directory with automatic method
using preset settings for *P vivax* with PvP01_v1 reference sequence, change to
base environment directory:

.. code-block:: console
  
      cd $NGSENV_BASEDIR

If running in an HPC/cluster system or workstation/server with 16-core or more,
use the following command to setup the base environment directory with full
version of PvP01_v1 setting:

.. code-block:: console

      bash <(curl -L https://raw.githubusercontent.com/vivaxgen/vgnpc-plasmodium-spp/main/Pvivax/PvP01_v1/setup.sh)

If running in a laptop or desktop with less than 16-core, use the following 
command to setup the base enviroment directory with lite version of PvP01_v1
setting:

.. code-block:: console

      bash <(curl -L https://raw.githubusercontent.com/vivaxgen/vgnpc-plasmodium-spp/main/Pvivax/PvP01_v1/setup-lite.sh)

The full version setup will take some time as it needs to download both the
PvP01_v1 genome sequence (~ 23MB), human GRCh38.p14 genome (~ 928MB),
uncompress the human genome, and generate index file for both PvP01 and the
human genome sequences using ``bwa-mem2``.

The lite version setup will take less time as it only needs to download the
PvP01_v1 genome sequence (~ 23MB) and generate index file for only PvP01
sequences.

.. note::

  The vivaxGEN github repository provides the list of available preset
  settings.
  However, if none of the preset settings are suitable, then the setup can be
  continued using manual method following steps described in
  :doc:`setup-base-env-dir`.


Running the Multi-Step Mode
---------------------------

This section of the tutorial shows the use of ``run-multistep-variant-caller``
single command, which provides the simple and quick way to perform multi-step
mode of the variant calling.
In this section, 2 samples of *P vivax* WGS data will be processed to get the
final result as a concatenated VCF file (a single VCF file containing all
chromosomes).

.. tip::

  For larger number of samples, it is advisable to have the final result as
  multiple VCF files, each contains a specific chromosome, since the downstream
  analysis then can be performed individually on each chromosome in parallel to
  speed up the analysis.

#.  Activate the environment by executing the ``activate`` script if the
    environment has not been activated::

	  /data/Pv-wgs/PvP01_v1/activate

#.  Enter the directory for containing data sets, and create a new directory,
    and enter to the new directory::

      cd $NGSENV_BASEDIR/sets
      mkdir my-tutorial
      cd my-tutorial

#.  Create a directory to hold the FASTQ read files::

	  mkdir reads-1

#.  Download read files related to 2 *P. vivax* sequence data from ENA (note
    that for working with public SRA read files, consider using
    `SRA-Repo <https://github.com/vivaxgen/sra-repo>`_ to manage and
    automatically download the read files)::

      cd reads-1
      wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR111/ERR111714/ERR111714_1.fastq.gz
      wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR111/ERR111714/ERR111714_2.fastq.gz
      wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR113/004/ERR1138854/ERR1138854_1.fastq.gz
      wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR113/004/ERR1138854/ERR1138854_2.fastq.gz
      cd ..

#.  Run the multi-step mode variant calling process by executing this single
    command::

      ngs-pl run-multistep-variant-caller -o batch-1 --paired reads-1/*.fastq.gz

    Wait until the process finishes.

#.  Inspect the ``batch-1`` directory by performing directorylisting::

      ls batch-1

    The following is the layout of the output directory:

    ``analysis/``
      This directory contains sample directory, eg. each sample and their
      associated files (input/output/log) are in their own directory.

    ``completed_samples/``
      This directory contains symbolic links to samples in ``analysis``
      directory that have been successfully called.
      The joint variant calling is performed only on samples in this
      directory.

    ``failed_samples/``
      This directory contains symbolic links to samples in ``analysis``
      that are failed during individual sample calling process.

    ``joint/``
      This directory contains all files pertinent to joint variant calling
      process.

    ``joint/concatenated.vcf.gz``
      This is the concatenated VCF file from chromosome-based VCF files
      inside ``joint/vcfs`` directory.
      This file is only available with ``--target concatenated_vcf`` option.

    ``joint/vcfs/``
      The final output of the joint variant calling is the per-chromosome
      VCF files in this directory.

    ``metafile/``
      This directory contains metafiles necessary for performing the whole
      variant calling process.
      Currently it holds the manifest file describing the sample name and its
      associated read files.

    ``reports/``
      This directory contains consolidated report files from completed samples
      in the ``completed_samples`` directory.
      Currently, it holds ``maps/`` directory (which links to BAM files of each
      samples) and ``depth-base/`` directory (which links to depth files
      generated by sambamba).

    ``stats.tsv``
      This file contains the statistics of each step of the process.

The main output file(s) of this whole variant calling process are VCF files
inside ``joint/vcfs`` and ``concatenated.vcf.gz``.

Now let assume that another batch of samples are available.
The following steps provide instructions to perform sample variant calling
and then do joint variant calling with the previous batch:

#.  Download read files related to another 2 of *P vivax* sequence data from
    SRA database::

      mkdir reads-2
      cd reads-2
      wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR527/ERR527357/ERR527357_1.fastq.gz
      wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR527/ERR527357/ERR527357_2.fastq.gz
      wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR152/ERR152414/ERR152414_1.fastq.gz
      wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR152/ERR152414/ERR152414_2.fastq.gz
      cd ..

#.  Run the multi-step variant calling with the new data, but only to the step
    of sample variant calling::

      ngs-pl run-multistep-variant-caller -o batch-2 --target sample_variant_calling reads-2/*.fastq.gz

    Wait until the process finishes.

#.  Run the joint-variant calling by combining the completed samples of
    ``batch-1`` and ``batch-2`` together::

      ngs-pl run-joint-variant-caller -o joint-batches --target concatenated_vcf batch-1/completed_samples batch-2/completed_samples

    Wait until the joint variant calling finishes.

#.  Inspect the directory ``joint-batches``.
    The per-chromosome VCF files would be in the ``joint-batches/vcfs``
    directory, while the concatenated VCF file containing all chromosomes in
    a single file would be ``joint-batches/concatenated.vcf.gz``.

Congratulation!
You now have sucessfully perform joint variant calling between 2 sample batches.


Working with SRA Data
---------------------

For working with many published FASTQ read files from SRA databases (NCBI SRA
or EMBL ENA), `SRA-Repo <https://github.com/vivaxgen/sra-repo>`_ can be used to
help downloading and managing SRA read files.

This part of tutorial requires ``SRA-Repo`` to be installed.
Follow the installation step in ``SRA-Repo`` github repository to install it
properly.

Open a new terminal/shell and change to the the tutorial directory.
Generate a tab-delimited sample file named ``my-samples.tsv`` with the content
as follow::

    SAMPLE      COUNTRY   SRA
    PH0098-C    C1        ERR216478,ERR490276
    PY0074-C    C2        ERR1138883

Activate SRA-Repo by activating its activation script, and fetch the SRA read
files in ``my-samples.tsv`` above::

    <YOUR_SRA_REPO_INSTALLATION>/bin/activate
    sra-repo.py fetch --ntasks 6 --samplefile my-samples.tsv:SAMPLE,SRA

The above command will download the SRA read files and store it inside the
``SRA-Repo`` installation directory.
After the download finishes, link the SRA read files to a new directory and
generate a manifest file::

    sra-repo.py link -o manifest-3.tsv --outdir reads-3 --samplefile my-samples.tsv:SAMPLE,SRA

In the terminal/shell with active NGS-Pipeline environment, perform sample
variant calling::

    ngs-pl run-multistep-variant-caller -o batch-3 --target sample_variant_calling -i manifest-3.tsv .

Note the dot (indicating current directory) at the last part of the above command.

Once the sample variant calling finishes, perform joint variant calling with the
previous batches::

    ngs-pl run-joint-variant-caller -o new-joint --target concatenated_vcf batch-1/completed_samples batch-2/completed_samples batch-3/completed_samples

Once the joint variant calling process finishes, inspect the result in the 
``new-joint``directory.


Expoloring Further
------------------

To read more about ``NGS-Pipeline`` features, please consult the rest of the
documentation.

