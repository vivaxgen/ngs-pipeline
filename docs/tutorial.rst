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
            refs/
            sets/

Note that only ``config.yaml`` is the mandatory name to be used for the config
file, while the rest of the file and directory names can be anything.
However, for consistency purposes, it is recommended to use the above file
and directory names.

The following are step-by-step instructions to in setting up the environment:

 1. Activate the vivaxGEN NGS-Pipeline environment by running its activation
    script, as noted after the automatic installation finished, eg::

      source NGS-PIPELINE_INSTALL_DIR/activate.sh

 2. Create the base working directory, eg: ``/data/Pvivax``::

      export BASEDIR=/data/Pvivax
	  mkdir $BASEDIR

 3. Generate activation script from the base directory::

      ngs-pl generate-activation-script -o $BASEDIR/activate.sh

 4. Edit the activation script as necessary, following the comments and notes
    in the activation script::

      vim $BASEDIR/activate.sh

 5. Prepare reference directory and populate the directory as necessary::

      mkdir $BASEDIR/refs


 6. Copy config file templates::

      ngs-pl generate-config-file -o $BASEDIR/config.yaml

 7. Edit the config file as necessary following the comments and notes in
    config file.


Running the Multi-Step Mode
---------------------------

* Activate the environment by sourcing the activation script::

	$ source BASEDIR/activate.sh
	


Running the Single-Step Mode
----------------------------




