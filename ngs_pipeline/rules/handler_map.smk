
# this file contains rules for processing mapped BAM files
# final output of this workflow is maps/mapped-final.bam,
# which is the merged BAM file of all input BAM files for that sample

# the mapping process in this snakemake file is:
# paired-maps -> filtered-maps -> deduped-maps (optional) -> merged-maps

from ngs_pipeline.rules import inc
from time import sleep

include: inc("reporter_stats.smk")

rule map_link:
    # this rule generate a hard link mapped-final.bam to the actual bam file
    # we use hard link because the actual bam might get deleted if set temporary
    localrule: True
    input:
        get_final_bam_files
    output:
        temp("<sp>maps/mapped-final-{idx}.bam")
    shell:
        "ln {input} {output}"



# to get only mapped, properly paired (FR) reads use filter-reads-orientation with
# --remove-unmapped --remove-RF --remove-FF --remove-RF --remove-trans

rule map_filter:
    threads: thread_allocations.get('map_filtering', 4)
    input:
        get_mapped_bam_file()
    output:
        temp("<sp>maps/mapped-filtered-{idx}.bam")
    log:
        log1 = "<sp>logs/filter_orientation-{idx}.log",
        log2 = "<sp>logs/samtools-sort-{idx}.log",
        read_orientation = "<sp>logs/read-orientation-{idx}.json"
    params:
        sample = get_sample,
        args = config.get('read_filters', '') or '--remove_unmapped',
    shell:
        "ngs-pl filter-reads-orientation --outstat {log.read_orientation} {params.args} {input} 2> {log.log1} "
        "| samtools sort -@4 -o {output} 2> {log.log2} "


rule map_dedup:
    threads: thread_allocations.get('map_deduping', 3)
    input:
        "<sp>maps/mapped-filtered-{idx}.bam"
    output:
        temp("<sp>maps/mapped-dedup-{idx}.bam")
    params:
        sample = get_sample,
    log:
        log1 = "<sp>logs/markdup-{idx}.log",
        markdup_stat = "<sp>logs/markdup-stat-{idx}.json",
    shell:
        "samtools markdup -r --json -f {log.markdup_stat} {input} {output} 2> {log.log1}"


rule map_merge_final:
    # this rule merges dedup input BAM
    threads: thread_allocations.get('map_merging', 4)
    input:
        #expand('<sp>maps/mapped-final-{idx}.bam', idx=IDXS)
        #get_merge_input_bam_files
        expand_sp('<sp>maps/mapped-final-{idx}.bam')
    output:
        temp_unless('<sp>maps/mapped-final.bam', keep_final_bam)
    params:
        sample = get_sample,
    run:
        if len(input) > 1:
            shell('samtools merge -@4 {output} {input}')
        else:
            # use hard link since input will be removed
            shell('ln {input} {output}')
        # to make index file newer by 1-sec to bam file
        sleep(2)
        shell('samtools index {output}')


rule collect_mapping_stats:
    threads: 1
    input:
        trims = expand_sp('<sp>logs/trimming_stat-{idx}.json'),
        maps = expand_sp('<sp>logs/{sample}-{idx}.stats.txt'),
        filtered = expand_sp('<sp>logs/mapped-filtered-{idx}.stats.txt'),
        dedups = expand_sp('<sp>logs/mapped-dedup-{idx}.stats.txt') if deduplicate else [],
        finals = expand_sp('<sp>logs/mapped-final-{idx}.stats.txt'),
        depths = expand_sp('<sp>logs/mapped-final-{idx}.depths.txt.gz'),
    params:
        sample = get_sample,
        trimmed = lambda wildcards, input: '--trimmed ' + ' --trimmed '.join(input.trims),
        mapped = lambda wildcards, input: '--mapped ' + ' --mapped '.join(input.maps),
        deduped = (lambda wildcards, input: '--dedup ' + ' --dedup '.join(input.dedups)) if deduplicate else '',
        finaled = lambda wildcards, input: '--final ' + ' --final '.join(input.finals),
        depthed = lambda wildcards, input: '--depth ' + ' --depth '.join(input.depths),
    output:
        '<sp>logs/stats.tsv'
    shell:
        'ngs-pl calculate-stats -o {output} --mindepth {min_depth} '
        '{params.trimmed} {params.mapped} {params.deduped} {params.finaled} {params.depthed} {params.sample}'

# EOF
