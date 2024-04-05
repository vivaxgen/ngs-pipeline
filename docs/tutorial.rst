TUTORIAL
========

The tutorial assumes that the pipeline has been installed as described in the
Quick Installation section of README.rst.
The tutorial will work through setting up a pipeline to process *Plasmodium
vivax* sequence data using PvP01_v1 reference sequence.


Base Enviroment directory
-------------------------

For normal workflows, vivaxGEN NGS-Pipeline requires that all necessary files,
with the exception of the pipeline and all supporting software themselves, to
reside in a single hierarchical directory (unless the pipeline has been set up
with specific settings).
The root of the hierarchical directory is called base environment directory,
and can be accessed with NGSENV_BASEDIR environment variabel.
All necessary files include the configuration, reference, data set (read files)
and the analysis result.
Since UNIX-based systems treat symbolic link files in identical ways as real
files, it is possible that the actual files reside in other directory and only
the symbolic links reside in the base environment directory.
This scheme was chosen to accomodate the cascading configuration mechanism,
where the pipeline process will walk across directory hierarchy to get
configuration file ``config.yaml`` from the working directory (or the sample
directory) to the base environment directory.

The recommended structure of the base environment directory is::

    NGSENV_BASEDIR/
                   activate
                   bashrc
                   config.yaml
                   configs/
                   refs/
                   sets/

Note that only ``bashrc``, ``config.yaml`` and ``configs`` are the mandatory
names to be used for the bash source file, main config file and configuration
directory, while the rest of the file and directory names can be anything.
However, for consistency purposes, it is recommended to use the above file
and directory names.
The ``refs`` directory is used to keep all reference files, including the
reference sequence and its index files, region files and any oher settings.
The ``sets`` directory is the main working directory to perform the analysis.


Preparing Base Environment Directory
------------------------------------

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

#.  Activate the environment by exectuing the ``activate`` script if the
    environment has not been activated::

	  /data/Pv-wgs/activate.sh

#.  Enter the directory for containing data sets, and create a new directory,
    and enter to the new directory::

      cd $NGSENV_BASEDIR/sets
      mkdir my-tutorial
      cd my-tutorial

#.  Create a directory to hold the FASTQ read files::

	  mkdir reads

#.  Download read files related to 2 *P. vivax* sequence data from ENA (note
    that for working with public SRA read files, consider using
    `SRA-Repo <https://github.com/vivaxgen/sra-repo>`_ to manage and
    automatically download the read files)::

      wget [url]


#.  Generate a manifest file from the read files::

      ngs-pl generate-manifest-file -o manifest.tsv reads/*.fastq.gz

#.  Edit the manifest file as necessary, such as changing the sample code.

#.  Run the sample directory preparation step (step-1)::

      ngs-pl prepare-sample-directory -o analysis -i manifest.tsv reads/

#.  Check for errors or warnings from the command.

#.  Run the sample variant-calling step (step-2)::

      ngs-pl run-sample-variant-caller analysis

#.  Run the joint variant-calling step (step-3)::

      ngs-pl run-joint-variant-caller -o joint analysis/

#.  Check the final VCF files in ``joint/vcfs`` directory.


Running the Single-Step Mode
----------------------------




