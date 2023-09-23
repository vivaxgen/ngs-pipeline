
# this is main snakemake module to perform sample variant calling producing
# individual GVCF file per sample

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
        'logs/stats.tsv',
        'logs/depths.png'


rule clean:
    shell:
        "rm -rf maps/ logs/ trimmed-reads/ gvcf/ .snakemake/*"


rule mapping:
    input:
        'maps/mapped-dedup.bam',
        'logs/stats.tsv'


include: config.get('reads_trimmer_wf', 'trimmer_fastp.smk')
include: config.get('reads_mapper_wf', 'mapper_bwa-mem2.smk')
include: config.get('base_calibrator_wf', 'calibratebase_gatk.smk')
include: config.get('variant_caller_wd', 'varcall_gatk.smk')

# the mapping process in this snakemake file is:
# paired-maps -> proper-maps -> deduped-maps -> merged-maps

# for wgs variant calling, we perform deduplication on mapped reads

rule map_proper:
    # this rule filter input BAM file for only mapped, properly paired (FR) reads
    threads: 3
    input:
        "maps/mapped-{idx}.bam"
    output:
        "maps/mapped-proper-{idx}.bam" if keep_proper_bam else temp("maps/mapped-proper-{idx}.bam")
    shell:
        "samtools view -F 0x4 -f 0x2 -q 15 -o  {output} {input}"

rule map_dedup:
    threads: 4
    input:
        "maps/mapped-proper-{idx}.bam"
    output:
        temp("maps/mapped-dedup-{idx}.bam")
    log:
        "logs/markdup-{idx}.log"
    shell:
        "samtools sort -@4 {input} | samtools markdup -r - {output} 2> {log}"


rule map_stats:
    threads: 1
    input:
        "maps/{filename}.bam"
    output:
        "logs/{filename}.stats.txt"
    shell:
        "samtools stats {input} > {output}"


rule map_merging_dedup:
    # this rule merges dedup input BAM
    threads: 8
    input:
        expand('maps/mapped-dedup-{idx}.bam', idx=IDXS)
    output:
        'maps/mapped-dedup.bam' if keep_deduplicated_bam else temp('maps/mapped-dedup.bam')
    run:
        if len(IDXS) > 1:
            shell('samtools merge -@8 {output} {input}')
        else:
            shell('cp {input} {output}')
        # to make index file newer by 1-sec to bam file
        sleep(2)
        shell('samtools index {output}')


rule map_depth:
    # this rule creates stats and depths of any .bam file using samtools stats
    threads: 1
    input:
        'maps/{filename}.bam'
    output:
        'logs/{filename}.depths.txt.gz'
    shell:
        'samtools depth {input} | gzip > {output}'


rule collect_stats:
    threads: 1
    input:
        trims = expand('logs/trimming_stat-{idx}.json', idx=IDXS),
        maps = expand('logs/mapped-{idx}.stats.txt', idx=IDXS),
        dedups = expand('logs/mapped-dedup-{idx}.stats.txt', idx=IDXS),
        depths = expand('logs/mapped-dedup-{idx}.depths.txt.gz', idx=IDXS),
    params:
        trimmed = lambda wildcards, input: '--trimmed ' + ' --trimmed '.join(input.trims),
        mapped = lambda wildcards, input: '--mapped ' + ' --mapped '.join(input.maps),
        deduped = lambda wildcards, input: '--dedup ' + ' --dedup '.join(input.dedups),
        depthed = lambda wildcards, input: '--depth ' + ' --depth '.join(input.depths),
    output:
        'logs/stats.tsv'
    shell:
        'collect_stats.py -o {output} {params.trimmed} {params.mapped} {params.deduped} {params.depthed} {sample}'


rule depth_plot:
    threads: 1
    input:
        "logs/mapped-dedup.depths.txt.gz"
    params:
        chroms = ('--chrom ' + ','.join(REGIONS)) if any(REGIONS) else ''
    output:
        'logs/depths.png'
    shell:
        'plot_depth.py --outplot {output} {params.chroms} --sort --infile {input} {sample}'


# EOF
