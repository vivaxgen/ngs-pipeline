TUTORIAL
========

The tutorial assumes that the pipeline has been installed as described in the
Quick Installation section of README.rst.
The tutorial will work through setting up a pipeline to process *Plasmodium
vivax* sequence data.

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
using preset settings for P vivax with PvP01_v1 reference sequence, change to
base environment directory and run the setup script::

      cd $NGSENV_BASEDIR
      bash <(curl PvP01_v1.sh)>


To continue peparing the base enviroment directory with manual method (eg, when
there is no preset settings from the vivaxGEN github repository), use the
following steps:

#.  Prepare the reference sequence in ``/data/Pv-wgs/refs`` directory.
    A P vivax reference sequence (PvP01) can be downloaded from PlasmoDB
    with the following commands::

      cd $NGSENV_BASEDIR/refs
      curl -O https://plasmodb.org/common/downloads/release-50/PvivaxP01/fasta/data/PlasmoDB-50_PvivaxP01_Genome.fasta
      ln -s PlasmoDB-50_PvivaxP01_Genome.fasta PvP01_v1.fasta

#.  Generate a YAML file denoting the sequence labels for the P vivax genomes,
    and edit the YAML output file accordingly to remove regions that are not
    to be analyzed (remove all regions started with ``Transfer``).
    The following command use ``sed`` to remove regions, but any text editor
    can be used instead as well::

      ngs-pl setup-references -o regions.yaml -f PvP01_v1.fasta
      sed -i.bak '/Transfer/d' regions.yaml

#.  Obtain a list of high-quality variants (chromosomes and base positions) of
    P vivax genome from other source.
    For P vivax, a bed file based on PASS variants of Pv4 data set from MalariaGEN
    project has been prepared and can be obtained using the following command:

    curl -O 
   

#.  NGS-Pipeline has the ability to filter out unwanted, or contaminant reads.
    As WGS of plasmodium parasite will also sequence the host (human) DNA,
    the host reference sequence can be downloaded as well (in this case,
    the human reference genome from NCBI ftp site)::

      cd $NGSENV_BASEDIR/refs
      curl -O https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
      gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz
      ln -s GCF_000001405.40_GRCh38.p14_genomic.fna GRCh38.p14.fasta

#.  Likewise, generate a YAML file for the human regions, but denote this
    as contaminants::

      ngs-pl setup-references -k contaminant_regions -o contaminants.yaml -f GRCh38.p14.fasta

#.  Concatenate both P vivax and human reference sequence to a single fasta file::

      cat PvP01_v1.fasta GRCh38.p14.fasta > PvP01_v1-GRCh38.p14.fasta

#.  Copy template config.yaml to the base environment directory and edit the config
    file as necessary, especially the path to the reference (the default values are
    suitable for many sequencing project)::

      cp $NGS_PIPELINE_BASE/config/config.yaml $NGENV_BASEDIR
      vim $NGSENV_BASEDIR/config.yaml

#.  Concatenate both ``regions.yaml`` and ``contaminants.yaml`` to config.yaml::

      cat regions.yaml contaminants.yaml >> $NGSENV_BASEDIR/config.yaml

#.  Check the configuration file::

      ngs-pl check-config-file $NGSENV_BASEDIR/config.yaml

    Fix any errors by editing the config.yaml, and then rerun the checking
    command until no more errors are reported.


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




