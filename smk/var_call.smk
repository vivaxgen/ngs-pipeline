
from time import sleep

# prepare necessary global parameters

include: "global_params.smk"

# prepare sample-related parameters

sample = pathlib.Path.cwd().name
IDXS, = glob_wildcards('reads/raw-{idx}_R1.fastq.gz')


# final output of this workflow

def get_final_file(w):
    return [f"gvcf/{sample}-{reg}.g.vcf.gz" for reg in REGIONS]


# define local rules
localrules: all, clean, mapping

rule all:
    input:
        get_final_file,
        'logs/mapped-dedup.stats.txt',
        'logs/stats.tsv'


rule clean:
    shell:
        "rm -rf maps/ logs/ trimmed-reads/ gvcf/ .snakemake/*"


rule mapping:
    input:
        'maps/mapped-dedup.bam',
        'logs/stats.tsv'


include: config.get('reads_trimmer_wf', 'trimmer_cutadapt.smk')
include: config.get('reads_mapper_wf', 'mapper_minimap2.smk')
include: config.get('base_calibrator_wf', 'calibratebase_gatk.smk')

# for wgs variant calling, we perform deduplication on mapped reads
rule map_dedup:
    threads: 4
    input:
        "maps/mapped-{idx}.bam"
    output:
        temp("maps/mapped-dedup-{idx}.bam")
    log:
        "logs/markdup-{idx}.log"
    shell:
        "samtools sort -@4 {input} | samtools markdup -r - {output} 2> {log}"


rule map_stats:
    threads: 1
    input:
        "maps/mapped-{idx}.bam"
    output:
        "logs/mapped-{idx}.stats.txt"
    shell:
        "samtools stats {input} > {output}"


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
            shell('cp {input} {output}')
        # to make index file newer by 1-sec to bam file
        sleep(2)
        shell('samtools index {output}')

include: "varcall_gatk.smk"


rule dedup_stats:
    threads: 1
    input:
        'maps/mapped-dedup.bam'
    output:
        'logs/mapped-dedup.stats.txt',
        'logs/mapped-dedup.depths.txt.gz'
    shell:
        'samtools stats {input} > {output[0]} && samtools depth {input} | gzip > {output[1]}'


rule stats:
    threads: 1
    input:
        trims = expand('logs/reads_trimming-{idx}.log', idx=IDXS),
        maps = expand('logs/mapped-{idx}.stats.txt', idx=IDXS),
        dedup = 'logs/mapped-dedup.stats.txt',
        depth = 'logs/mapped-dedup.depths.txt.gz'
    params:
        trimmed = lambda wildcards, input: '--trimmed ' + ' --trimmed '.join(input.trims),
        mapped = lambda wildcards, input: '--mapped ' + ' --mapped '.join(input.maps)
    output:
        'logs/stats.tsv'
    shell:
        'collect_stats.py -o {output} {params.trimmed} {params.mapped} --dedup {input.dedup} --depth {input.depth} {sample}'


# EOF
