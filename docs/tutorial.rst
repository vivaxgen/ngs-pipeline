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



Running the Multi-Step Mode
---------------------------

* Activate the environment by sourcing the activation script::

	$ source BASEDIR/activate.sh
	


Running the Single-Step Mode
----------------------------




