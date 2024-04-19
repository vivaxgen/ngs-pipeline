
Setting Up Base Enviroment Directory
====================================

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


Manual Setting Up of Base Environment Directory
-----------------------------------------------

The base environment directory can be setup automatically using preset
configurations available from vivaxGEN github repository.
However, if none of the preset configuration is suitable for processing the
data (such as the species or reference sequences are not in the list of preset
configuration), then a manual method to prepare the base environment directory
is necessary (which is what this document is about).

Before setting up the environment, decide on which reference sequence to be
used and where to obtain it.
Optionally, a list of known variants based on that reference sequence can be
used to calibrate the base calling of the data.
Even though NGS-Pipeline can still work without the list of known variants,
it is recommended to obtain the variant list from other documents or sources.

The following are the step-by-step instruction on how to prepare the setting of
base environment directory using *P vivax* PvP01 reference sequence as an
example.

Perform the initial step to create the base environment directory as follows.
The ``NGS-PIPELINE_INSTALL_DIR`` is assumed to be the directory where the
NGS-Pipeline is installed, while ``/data/Pv-wgs`` is the new base environment
directory::

      NGS-PIPELINE_INSTALL_DIR/bin/activate
      ngs-pl setup-base-directory /data/Pv-wgs
      exit
      # edit /data/Pv-wgs/activate if necessary
      /data/Pv-wgs/activate

Once the new base environment has been activated, NGSENV_BASEDIR environment
variable will be accessible.
The next steps are as follow:

#.  Create the directories to hold the reference sequence and all related files
    and change to the directory:

      mkdir -p $NGSENV_BASEDIR/refs/PvP01_v1/known_variants
      cd $NGSENV_BASEDIR/refs/PvP01_v1/

#.  Prepare the reference sequence in ``/data/Pv-wgs/refs/PvP01_v1`` directory.
    A *P vivax* reference sequence (PvP01_v1) can be downloaded from PlasmoDB
    with the following commands::

      curl -O https://plasmodb.org/common/downloads/release-50/PvivaxP01/fasta/data/PlasmoDB-50_PvivaxP01_Genome.fasta
      ln -s PlasmoDB-50_PvivaxP01_Genome.fasta PvP01_v1.fasta

#.  Generate a YAML file denoting the sequence labels for the ^P vivax* genomes,
    and edit the YAML output file accordingly to remove regions that are not
    to be analyzed (remove all regions started with ``Transfer``).
    The following command use ``sed`` to remove regions, but any text editor
    can be used instead as well::

      ngs-pl setup-references -o regions.yaml -f PvP01_v1.fasta
      sed -i.bak '/Transfer/d' regions.yaml

#.  Obtain a list of high-quality variants (chromosomes and base positions) of
    *P vivax* genome from other source.
    For *P vivax*, a bed file based on PASS variants of Pv4 data set from
    MalariaGEN project has been prepared in BED format and can be obtained
    using the following command:

	curl -O https://raw.githubusercontent.com/vivaxgen/vgnpc-plasmodium-spp/main/Pvivax/PvP01_v1/known-variants.bed.gz

#.  Split the variants based on their chromosome names into disctint files
    to speed up the base calibration process:

      python3 -c "import yaml; [open(f'known-variants/{reg}.bed', 'w') for reg in yaml.safe_load(open('regions.yaml'))['regions']]"
      zcat known-variants.bed.gz | awk '{print > "known-variants/" $1 ".bed"}'
	for fn in known-variants/*.bed; do echo "Processing ${fn}"; bgzip -f ${fn}; tabix ${fn}.gz; done

#.  NGS-Pipeline has the ability to filter out unwanted, or contaminant reads.
    As WGS of plasmodium parasite will also sequence the host (human) DNA,
    the host reference sequence can be downloaded as well (in this case,
    the human reference genome from NCBI ftp site)::

      curl -O https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
      gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz
      ln -s GCF_000001405.40_GRCh38.p14_genomic.fna GRCh38.p14.fasta

#.  Likewise, generate a YAML file for the human regions, but denote this
    as contaminants::

      ngs-pl setup-references -n -k contaminant_regions -o contaminants.yaml -f GRCh38.p14.fasta

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
