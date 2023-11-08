
# prepare reference 
import os
import pathlib
from time import sleep

ngs_pipeline_basedir = os.environ['NGS_PIPELINE_BASE']
ngsenv_basedir = os.environ['NGSENV_BASEDIR']

refmap = ngsenv_basedir + '/' + config.get('refmap_file', 'NOFILE')
refseq = ngsenv_basedir + '/' + config['refseq_file']

# prepare sample-related parameters

sample = pathlib.Path.cwd().name
IDXS, = glob_wildcards('reads/raw-{idx}_R0.fastq.gz')


rule all:
    localrule: True
    input:
        'maps/mapped-final.bam'


rule clean:
    localrule: True
    shell:
        "rm -rf maps/ logs/ trimmed-reads/ gvcf/ .snakemake/*"


rule trim_reads:
    threads: 4
    input:
        'reads/raw-{idx}_R0.fastq.gz'
    output:
        'trimmed-reads/trimmed-{idx}_R0.fastq.gz'
    shell:
        "gunzip -c {input} | chopper -q 10 -l 1000 --headcrop 30 --tailcrop 30 | gzip > {output}"


rule mapping:
    threads: 8
    input:
        'trimmed-reads/trimmed-{idx}_R0.fastq.gz'
    output:
        temp('maps/mapped-sorted-{idx}.bam')
    params:
        rg = lambda w: f"-R @RG\\\\tID:{sample}\\\\tSM:{sample}\\\\tLB:LIB-{sample}",
    shell:
        "minimap2 -a {refmap} {params.rg} {input} "
        "|  samtools sort -@4 -o {output} "


rule map_merge_final:
    # this rule merges dedup input BAM
    threads: 4
    localrule: True
    input:
        expand('maps/mapped-sorted-{idx}.bam', idx=IDXS)
    output:
        'maps/mapped-final.bam'
    run:
        if len(IDXS) > 1:
            shell('samtools merge -@4 {output} {input}')
        else:
            # use hard link since input will be removed
            shell('ln {input} {output}')
        # to make index file newer by 1-sec to bam file
        sleep(2)
        shell('samtools index {output}')


# EOF
