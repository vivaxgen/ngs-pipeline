
# this is main snakemake module to perform sample variant calling producing
# BAM from ONT long reads

from time import sleep

# prepare necessary global parameters

include: "global_params.smk"

include: "utilities.smk"
# prepare sample-related parameters

sample = pathlib.Path.cwd().name
IDXS, = glob_wildcards('reads/raw-{idx}_R0.fastq.gz')

interval_file = config.get('interval_file', None)
interval_dir = config.get('interval_dir', None)

def get_interval(w):
    global interval_dir, interval_file
    if interval_dir:
        return f'--targets {interval_dir}/{w.reg}.bed'
    if interval_file:
        return f'--targets {interval_file}'
    return f'--region {w.reg}'

# utilities

def get_final_bam(w):
    """ return the final bam file for further processing """
    #return ("maps/mapped-dedup-{idx}.bam" if deduplicate else "maps/mapped-filtered-{idx}.bam")
    return ("maps/mapped-filtered-{idx}.bam")


def get_final_file(w):
    """ return final file """
    return "maps/mapped-final.bam"


include: config.get('reads_trimmer_wf', 'ssf_trimmer_fastplong.smka')
include: config.get('reads_mapper_wf', 'ssf_mapper_minimap2_lr.smk')
#include: config.get('variant_caller_wf', 'ssf_varcall_gatk.smk')
include: config.get('stats_wf', 'ssf_stats.smk')

# final output of this workflow


rule all:
    localrule: True
    input:
        get_final_file,
        'logs/mapped-final.stats.txt',
        'logs/mapped-final.depth-base.tsv.gz',
        'logs/stats.tsv',
        'logs/depths.png'


rule mapping:
    localrule: True
    input:
        'maps/mapped-final.bam',
        'logs/mapped-final.stats.txt',
        'logs/mapped-final.depth-base.tsv.gz',
        'logs/stats.tsv',
        #'logs/depths.png',

rule clean:
    localrule: True
    shell:
        "rm -rf maps/ logs/ trimmed-reads/ gvcf/ .snakemake/*"


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


assert deduplicate is False, "deduplication not supported for ONT reads"

rule collect_stats:
    threads: 1
    input:
        trims = expand('logs/trimming_stat-{idx}.json', idx=IDXS),
        maps = expand('logs/mapped-sorted-{idx}.stats.txt', idx=IDXS),
        filtered = expand('logs/mapped-sorted-{idx}.stats.txt', idx=IDXS),
        dedups = expand('logs/mapped-dedup-{idx}.stats.txt', idx=IDXS) if deduplicate else [],
        finals = expand('logs/mapped-sorted-{idx}.stats.txt', idx=IDXS),
        depths = expand('logs/mapped-sorted-{idx}.depths.txt.gz', idx=IDXS),
    params:
        sample = sample,
        trimmed = lambda wildcards, input: '--trimmed ' + ' --trimmed '.join(input.trims),
        mapped = lambda wildcards, input: '--mapped ' + ' --mapped '.join(input.maps),
        deduped = (lambda wildcards, input: '--dedup ' + ' --dedup '.join(input.dedups)) if deduplicate else '',
        finaled = lambda wildcards, input: '--final ' + ' --final '.join(input.finals),
        depthed = lambda wildcards, input: '--depth ' + ' --depth '.join(input.depths),
    output:
        'logs/stats.tsv'
    shell:
        'ngs-pl calculate-stats -o {output} --mindepth {min_depth} '
        '{params.trimmed} {params.mapped} {params.deduped} {params.finaled} {params.depthed} {sample}'


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
        fb_args = "--haplotype-length -1 --min-coverage 10 " + config["freebayes_extra_flags"]
    shell:
        "freebayes --fasta-reference {refseq} --ploidy {ploidy} --min-alternate-count 2 {params.reg} "
        "{params.fb_args} {input} 2> {log.fb} "
        "| bcftools view -o {output} "


# EOF
