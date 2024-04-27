vivaxGEN NGS-Pipeline Documentation
===================================

vivaxGEN NGS-Pipeline is an open-source, non-opiniated pipeline for variant
calling (upstream/secondary processing) for paired-end short reads or singleton
long reads NGS data.

As an non-opiniated pipeline, vivaxgen NGS-Pipeline supports the following
program/system:

* read trimmers: fastp, cutadapt, chopper

* mappers: bwa/bwa-mem2, bowtie2, minimap2

* variant callers: GATK4 (following GATK Best Practises as much as possible),
  FreeBayes, BCFTools mpileup/call, Clair3

vivaxGEN NGS-Pipeline can be used for variant discovery and targeted/panel
variant calling.
Variant discovery is suitable for various studies such as population genetics
and GWAS, usually from WGS data.
Targeted and panel variant calling is suitable for amplicon-sequencing data
and experiments that do not require joint-variant calling such as single sample
analysis for detecting status of certain variants (reference/alternate alleles,
heterozygotes or no data).

vivaxGEN NGS-Pipeline can be installed on laptops, servers or cluster/HPC
system without the need of administrator/root privileges.
The only requirement is a UNIX-based system supported by `micromamba
<https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html>`
(eg. various Linux distributions, WSL2, MacOSX) with preinstalled ``curl``
and ``bash`` (these two programs are usually installed in the base system).

