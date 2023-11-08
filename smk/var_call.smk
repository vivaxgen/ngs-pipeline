
# this is main snakemake module to perform sample variant calling producing
# individual GVCF file per sample

from time import sleep

# prepare necessary global parameters

include: "global_params.smk"

# prepare sample-related parameters

sample = pathlib.Path.cwd().name
IDXS, = glob_wildcards('reads/raw-{idx}_R1.fastq.gz')


# utilities

def get_final_bam(w):
    """ return the final bam file for further processing """
    return ("maps/mapped-dedup-{idx}.bam" if deduplicate else "maps/mapped-filtered-{idx}.bam")


# final output of this workflow

def get_final_file(w):
    return [f"gvcf/{sample}-{reg}.g.vcf.gz" for reg in REGIONS]


# define local rules
localrules: all, clean, mapping


rule all:
    input:
        get_final_file,
        'logs/mapped-final.stats.txt',
        'logs/stats.tsv',
        'logs/depths.png'


rule clean:
    shell:
        "rm -rf maps/ logs/ trimmed-reads/ gvcf/ .snakemake/*"


rule mapping:
    input:
        'maps/mapped-final.bam',
        'logs/stats.tsv'


include: config.get('reads_trimmer_wf', 'trimmer_fastp.smk')
include: config.get('reads_mapper_wf', 'mapper_bwa-mem2.smk')
include: config.get('base_calibrator_wf', 'calibratebase_gatk.smk')
include: config.get('variant_caller_wf', 'varcall_gatk.smk')

# the mapping process in this snakemake file is:
# paired-maps -> proper-maps -> deduped-maps -> merged-maps

# for wgs variant calling, we perform deduplication on mapped reads

rule map_link:
    # this rule generate a hard link mapped-final.bam to the actual bam file
    # we use hard link because the actual bam might get deleted if set temporary
    localrule: True
    input:
        get_final_bam
    output:
        temp("maps/mapped-final-{idx}.bam")
    shell:
        "ln {input} {output}"


rule map_proper:
    # this rule filter input BAM file for only mapped, properly paired (FR) reads
    threads: 3
    input:
        "maps/mapped-{idx}.bam"
    output:
        "maps/mapped-proper-{idx}.bam" if keep_proper_bam else temp("maps/mapped-proper-{idx}.bam")
    shell:
        "samtools view -F 0x4 -f 0x2 -q 15 -o  {output} {input}"


rule map_filter:
    threads: 4
    input:
        "maps/mapped-{idx}.bam"
    output:
        temp("maps/mapped-filtered-{idx}.bam")
    log:
        log1 = "logs/filter_orientation-{idx}.log",
        log2 = "logs/samtools-sort-{idx}.log",
        read_orientation = "logs/read-orientation-{idx}.json"
    params:
        args = config.get('read_filters', '') or '--remove_unmapped'
    shell:
        "filter_reads_orientation.py --outstat {log.read_orientation} {params.args} {input} 2> {log.log1} "
        "| samtools sort -@4 -o {output} 2> {log.log2} "


rule map_dedup:
    threads: 3
    input:
        "maps/mapped-filtered-{idx}.bam"
    output:
        temp("maps/mapped-dedup-{idx}.bam")
    log:
        log1 = "logs/markdup-{idx}.log",
        markdup_stat = "logs/markdup-stat-{idx}.json",
    shell:
        "samtools markdup -r --json -f {log.markdup_stat} {input} {output} 2> {log.log1}"


rule map_stats:
    threads: 1
    input:
        "maps/{filename}.bam"
    output:
        "logs/{filename}.stats.txt"
    shell:
        "samtools stats {input} > {output}"


rule map_merge_final:
    # this rule merges dedup input BAM
    threads: 4
    input:
        expand('maps/mapped-final-{idx}.bam', idx=IDXS)
    output:
        'maps/mapped-final.bam' if keep_final_bam else temp('maps/mapped-final.bam')
    run:
        if len(IDXS) > 1:
            shell('samtools merge -@4 {output} {input}')
        else:
            # use hard link since input will be removed
            shell('ln {input} {output}')
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
        filtered = expand('logs/mapped-filtered-{idx}.stats.txt', idx=IDXS),
        dedups = expand('logs/mapped-dedup-{idx}.stats.txt', idx=IDXS) if deduplicate else [],
        finals = expand('logs/mapped-final-{idx}.stats.txt', idx=IDXS),
        depths = expand('logs/mapped-final-{idx}.depths.txt.gz', idx=IDXS),
    params:
        trimmed = lambda wildcards, input: '--trimmed ' + ' --trimmed '.join(input.trims),
        mapped = lambda wildcards, input: '--mapped ' + ' --mapped '.join(input.maps),
        deduped = (lambda wildcards, input: '--dedup ' + ' --dedup '.join(input.dedups)) if deduplicate else '',
        finaled = lambda wildcards, input: '--final ' + ' --final '.join(input.finals),
        depthed = lambda wildcards, input: '--depth ' + ' --depth '.join(input.depths),
    output:
        'logs/stats.tsv'
    shell:
        'collect_stats.py -o {output} --mindepth {min_depth} '
        '{params.trimmed} {params.mapped} {params.deduped} {params.finaled} {params.depthed} {sample}'


rule depth_plot:
    threads: 1
    input:
        "logs/mapped-final.depths.txt.gz"
    params:
        chroms = ('--chrom ' + ','.join(REGIONS)) if any(REGIONS) else ''
    output:
        'logs/depths.png'
    shell:
        'plot_depth.py --outplot {output} {params.chroms} --sort --infile {input} {sample}'


# EOF
