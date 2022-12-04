
# prepare necessary parameters

import os

ngs_pipeline_basedir = os.environ['NGS_PIPELINE_BASE']
ngsenv_basedir = os.environ['NGSENV_BASEDIR']

refmap = ngsenv_basedir + '/' + config.get('refmap_file', 'NOFILE')
refseq = ngsenv_basedir + '/' + config['refseq_file']
knownsites_file = ngsenv_basedir + '/' + config.get('knownsites_file', '')

java_opts = ''

# check available read files

wildcard_constraints:
    idx = '\\d+'

IDXS, = glob_wildcards('reads/raw-{idx}_R1.fastq.gz')

# prepare global params

include: "global_params.smk"

# final output of this workflow

def get_final_file(w):
    return [ f"gvcf/{sample}-{reg}.g.vcf" for reg in REGIONS ]

rule all:
    input:
        get_final_file


rule clean:
    shell:
        "rm -rf maps/ logs/ trimmed-reads/ gvcf/ .snakemake"


include: config.get('reads_trimmer_wf', 'trimmer_cutadapt.smk')
include: config.get('reads_mapper_wf', 'mapper_minimap2.smk')
include: config.get('base_calibrator_wf', 'calibratebase_gatk.smk')

# for wgs variant calling, we perform deduplication on mapped reads
rule map_dedup:
    threads: 4
    input:
        "maps/mapped-{idx}.bam"
    output:
        "maps/mapped-dedup-{idx}.bam"
    shell:
        "samtools sort -@4 {input} | samtools markdup -r - {output}"

rule map_merging:
    threads: 8
    input:
        expand('maps/mapped-dedup-{idx}.bam', idx=IDXS)
    output:
        'maps/mapped-dedup.bam'
    run:
        if len(IDXS) > 1:
            shell('samtools merge -@8 {output} {input}')
        else:
            shell('ln -sr {input} {output}')
        shell('samtools index {output}')

include: "varcall_gatk.smk"

# EOF
