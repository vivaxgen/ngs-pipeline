# msf_merge_map.smk - ngs-pipeline rules
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2023, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

# required params
# - read_files


rule map_filter_orientation:
    # this rule filter the mapped reads based on orientation and then sorted by coordinate
    threads: thread_allocations.get('map_filtering', 4)
    input:
        bam = "{pfx}/{sample}/maps/{sample}-{idx}.bam",
    output:
        bam = temp("{pfx}/{sample}/maps/mapped-filtered-{idx}.bam")
    log:
        log1 = "{pfx}/{sample}/logs/filter_orientation-{idx}.log",
        log2 = "{pfx}/{sample}/logs/samtools-sort-{idx}.log",
        read_orientation = "{pfx}/{sample}/logs/read-orientation-{idx}.json"
    params:
        args = config.get('read_filters', '') or '--remove_unmapped',
    shell:
        "ngs-pl filter-reads-orientation --outstat {log.read_orientation} {params.args} {input} 2> {log.log1} "
        "| samtools sort -@4 -o {output} 2> {log.log2} "


def get_sorted_bam_files(w):
    return expand('{{pfx}}/{{sample}}/maps/mapped-filtered-{idx}.bam', idx=read_files.get_indexes(w.sample))

rule msf_merge_map:
    # this rule merges all bams into a single sorted bam files
    # note: stamtools merge sort the positions unless invoked with -n
    threads: 4
    input:
        #expand('{{pfx}}/{{sample}}/maps/sorted-{idx}.bam', idx=read_files.get_indexes)
        get_sorted_bam_files
    output:
        bam = "{pfx}/{sample}/maps/sorted.bam"
    run:
        if len(input) > 1:
            shell('samtools merge -@4 {output.bam} {input}')
        else:
            # use hard link since input will be removed
            shell('ln {input} {output}')


rule msf_filter_target:
    # this rule filter the mapped reads based on target regions
    threads: thread_allocations.get('map_filter_target', 4)
    input:
        bam = "{pfx}/{sample}/maps/sorted.bam",
    output:
        bam = temp("{pfx}/{sample}/maps/target.bam")
    params:
        region_opts = f'-L {targetregion_file}' if targetregion_file else ""
    shell:
        "samtools view -o {output.bam} {params.region_opts} {input.bam}"


# EOF
