
# prepare reference 
import os
import pathlib
from time import sleep

ngs_pipeline_basedir = os.environ['NGS_PIPELINE_BASE']
ngsenv_basedir = os.environ['NGSENV_BASEDIR']

refmap = ngsenv_basedir + '/' + config.get('refmap_file', 'NOFILE')
refseq = ngsenv_basedir + '/' + config['refseq_file']

ploidy = config.get('ploidy', 2)
freebayes_extra_flags = config.get('freebayes_extra_flags', '')

PARTIALS = None
REGIONS = config.get('regions', [])
if isinstance(REGIONS, dict):
    PARTIALS = REGIONS
    REGIONS = PARTIALS.keys()

interval_file = config.get('interval_file', None)
interval_dir = config.get('interval_dir', None)

def get_interval(w):
    global interval_dir, interval_file
    if interval_dir:
        return f'--targets {interval_dir}/{w.reg}.bed'
    if interval_file:
        return f'--targets {interval_file}'
    return f'--region {w.reg}'

# prepare sample-related parameters

sample = pathlib.Path.cwd().name
IDXS, = glob_wildcards('reads/raw-{idx}_R0.fastq.gz')

include: "utilities.smk"


# final output of this workflow

def get_final_file(w):
    return [f"gvcf/{sample}-{reg}.g.vcf.gz.csi" for reg in REGIONS]


rule all:
    localrule: True
    input:
        get_final_file

rule mapping:
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


rule map_reads:
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


rule fb_varcall:
    threads: 1
    input:
        'maps/mapped-final.bam'
    output:
        f'gvcf/{sample}-{{reg}}.g.vcf.gz'
    log:
        fb = 'logs/freebayes-{reg}.log'
    params:
        reg = get_interval,
        fb_args = "--haplotype-length -1 --min-coverage 10 " + freebayes_extra_flags
    shell:
        "freebayes --fasta-reference {refseq} --ploidy {ploidy} --min-alternate-count 2 {params.reg} "
        "{params.fb_args} {input} 2> {log.fb} "
        "| bcftools view -o {output} "
# EOF
