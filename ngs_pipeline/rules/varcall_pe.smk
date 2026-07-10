# varcall_pe.smk
# variant calling workflow for paired-end data

from ngs_pipeline import cerr
from ngs_pipeline.rules import inc
from time import sleep

# utilities


# define local rules
localrules: all, clean, mapping

include: inc("handler_map.smk")

rule all:
    input:
        get_final_file,
        add_prefix('logs/mapped-final.stats.txt'),
        add_prefix('logs/mapped-final.depth-base.tsv.gz'),
        add_prefix('logs/stats.tsv'),
        add_prefix('logs/depths.png'),

rule all_no_qc:
    input:
        get_final_file,

rule clean:
    shell:
        "rm -rf maps/ logs/ trimmed-reads/ gvcf/ .snakemake/*"


rule trimming:
    input:
        #expand(add_prefix('trimmed-reads/trimmed-{idx}_R1.fastq.gz'), idx=IDXS),
        #expand(add_prefix('trimmed-reads/trimmed-{idx}_R2.fastq.gz'), idx=IDXS),
        expand_sp('<sp>trimmed-reads/trimmed-{idx}_R1.fastq.gz'),
        expand_sp('<sp>trimmed-reads/trimmed-{idx}_R2.fastq.gz')


rule mapping:
    input:
        expand_sp('maps/mapped-final.bam'),
        add_prefix('logs/mapped-final.stats.txt'),
        add_prefix('logs/mapped-final.depth-base.tsv.gz'),
        add_prefix('logs/stats.tsv'),
        add_prefix('logs/depths.png'),


# default is fastp -> bwa-mem2 -> GATK

#include: config.get('reads_trimmer_wf', 'trimmer_fastp.smk')
#include: config.get('reads_mapper_wf', 'mapper_bwa-mem2.smk')
#include: config.get('variant_caller_wf', 'varcall_gatk.smk')
#include: config.get('stats_wf', 'stats.smk')

include: inc(config.get("trimmer_wf", "trimmer_fastp.smk"))

include: inc(config.get("mapper_wf", "mapper_bwa-mem2.smk"))

include: inc(config.get("varcall_wf", "varcall_gatk.smk"))

# EOF
