

rule map_filter:
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

# EOF