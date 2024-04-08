TUTORIAL
========

The tutorial assumes that the pipeline has been installed as described in the
Quick Installation section of README.rst.
The tutorial will work through setting up a pipeline to process *Plasmodium
vivax* sequence data using PvP01_v1 reference sequence.


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
    script, as noted after the automatic installation finished, eg::

      NGS-PIPELINE_INSTALL_DIR/activate

#.  Setup the base working directory, eg: ``/data/Pv-wgs``::

      ngs-pl setup-base-directory /data/Pv-wgs

#.  Change to ``/data/Pv-wgs`` directory and edit the ``activate`` script as
    necessary (eg, changing the prompt notification).
    Note that any text editor (eg: nano, vim, etc) can be used::

      cd /data/Pv-wgs
      vi activate

#.  Exit and activate the environment using the new ``activate`` file::

      exit
      /data/Pv-wgs/activate

    Once activated, the environment directory can be accessed using environment
    variable ``NGSENV_BASEDIR``.

To continue preparing the base enviroment directory with automatic method
using preset settings for *P vivax* with PvP01_v1 reference sequence, change to
base environment directory and run the setup script::

      cd $NGSENV_BASEDIR
      bash <(curl -L https://raw.githubusercontent.com/vivaxgen/vgnpc-plasmodium-spp/main/Pvivax/PvP01_v1/setup.sh)>

The above step will take some time as it needs to download both the PvP01 genome
sequence (~ 23MB), human GRCh38.p14 genome (~ 928MB), uncompress the human genome,
and generate index file for both PvP01 and the human genome sequences.

The vivaxGEN github repository provides the list of available preset settings.
However, if none of the preset settings are suitable, then the setup can be
continued using manual method following steps described in
`this document <setup-base-env-dir.rst>`_.

Running the Multi-Step Mode
---------------------------

#.  Activate the environment by executing the ``activate`` script if the
    environment has not been activated::

	  /data/Pv-wgs/activate.sh

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
      wget [url]
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

    ``concatenated.vcf.gz``
      This is the concatenated VCF file from chromosome-based VCF files
      inside ``joint/vcfs`` directory.

    ``failed_samples/``
      This directory contains symbolic links to samples in ``analysis``
      that are failed during individual sample calling process.

    ``joint/``
      This directory contains all files related to joint variant calling
      process.
      The final output of the joint variant calling is the per-chromosome
      VCF files in ``joint/vcfs/`` directory, which is being concatenated
      as ``concatenated.vcf.gz``.

    ``metafile``
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
and then do joint variant calling with the previous batch::

#.  Download read files related to another 2 of *P vivax* sequence data from
    SRA database::

      mkdir reads-2
      cd reads-2
      wget [urls]
      cd ..

#.  Run the multi-step variant calling with the new data, but only to the step
    of sample variant calling::

      ngs-pl run-multistep-variant-caller -o batch-2 --target GVCF reads-2/*.fastq.gz

    Wait until the process finishes.

#.  Run the joint-variant calling by combining the completed samples of
    ``batch-1`` and ``batch-2`` together::

      ngs-pl run-joint-variant-caller -o joint-batches --target concatenated_vcf batch-1/completed_samples batch-2/completed_samples

    Wait until the joint variant calling finishes.

#.  Inspect the directory ``joint-batches``.
    The per-chromosome VCF files would be in the ``joint-batches/vcfs``
    directory, while the concatenated VCF file containing all chromosomes in
    a single file would be ``joint-batches/concatenated.vcf.gz``.


Working with SRA Data
---------------------

For working with many published FASTQ read files from SRA databases (NCBI SRA
or EMBL ENA), `SRA-Repo <https://github.com/vivaxgen/sra-repo>` can be used to
help downloading and managing 