TUTORIAL
========

The tutorial assumes that the pipeline has been installed as described in the
Quick Installation section of README.rst. The tutorial will work through
setting up a pipeline to process *Plasmodium vivax* sequence data.


Setting-up Environment
----------------------

To set up a working environment, vivaxGEN ngs-pipeline requires that all
necessary files reside in a single hierarchical directory.
Since UNIX-based systems treat symbolic link files in identical ways as real
files, it is possible that the actual files reside in other directory and only
the symbolic links reside in the vivaxGEN ngs-pipeline environment directory.

The recommended structure of the base working directory is::

    BASEDIR/
            activate.sh
            config.yaml
            configs/
            refs/
            sets/

Note that only ``config.yaml`` and ``configs`` are the mandatory name to be
used for the config file and directory, while the rest of the file and
directory names can be anything.
However, for consistency purposes, it is recommended to use the above file
and directory names.

The following are step-by-step instructions to in setting up the environment:

1.  Activate the vivaxGEN NGS-Pipeline environment by running its activation
    script, as noted after the automatic installation finished, eg::

      source NGS-PIPELINE_INSTALL_DIR/activate.sh

2.  Create the base working directory, eg: ``/data/Pvivax``::

      export BASEDIR=/data/Pvivax
      mkdir $BASEDIR

3.  Generate activation script from the base directory::

      ngs-pl generate-activation-script -o $BASEDIR/activate.sh

4.  Edit the activation script as necessary, following the comments and notes
    in the activation script::

      vim $BASEDIR/activate.sh

5.  Activate the new activation script::

      source $BASEDIR/activate.sh

6.  Prepare reference directory and populate the directory as necessary::

      mkdir $BASEDIR/refs

7.  Copy config file templates::

      ngs-pl generate-config-file -o $BASEDIR/config.yaml

8.  Edit the config file as necessary following the comments and notes in
    config file::

      vim $BASEDIR/config.yaml

9.  Check the configuration file::

      ngs-pl check-config-file $BASEDIR/config.yaml

    Fix any errors by editing the config.yaml, and then rerun the checking
    command until no more errors are reported.


Running the Multi-Step Mode
---------------------------

#.  Open a new shell, and activate the environment by sourcing the activation
    script::

	  source /data/Pvivax/activate.sh

#.  Enter the directory for containing data sets, and create a new directory,
    and enter to the new directory::

      cd /data/Pvivax/sets
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




