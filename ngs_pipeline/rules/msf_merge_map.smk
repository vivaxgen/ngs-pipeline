# msf_merge_map.smk - ngs-pipeline rules
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2023, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

# required params
# - read_files

sam_flt = config.get('sam_read_filters', '')

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
        additional_sam_filter = f"-e \"{sam_flt}\"" if sam_flt != '' else "",
    shell:
        "ngs-pl filter-reads-orientation --outstat {log.read_orientation} {params.args} {input} 2> {log.log1} "
        "| samtools sort -@4 2> {log.log2} "
        "| samtools view {params.additional_sam_filter} -o {output} "


rule msf_final_map:
    # this rule prepares final bam file by filtering the mapped reads based on target regions
    threads: thread_allocations.get('map_filter_target', 4)
    input:
        bam = "{pfx}/{sample}/maps/mapped-filtered-{idx}.bam",
    output:
        bam = temp("{pfx}/{sample}/maps/mapped-final-{idx}.bam")
    params:
        region_opts = f'-L {targetregion_file}' if targetregion_file else ""
    run:
        # if no target region specified, just symbolic link the input as output
        if not params.region_opts:
            shell(f"ln -srf {input.bam} {output.bam}")
        else:
            shell(f"samtools view -o {output.bam} {params.region_opts} {input.bam}")


def get_final_bam_files(w):
    return expand('{{pfx}}/{{sample}}/maps/mapped-final-{idx}.bam', idx=read_files.get_indexes(w.sample))


rule msf_merge_map:
    # this rule merges all bams into a single sorted bam files
    # note: stamtools merge sort the positions unless invoked with -n
    threads: 4
    input:
        #expand('{{pfx}}/{{sample}}/maps/sorted-{idx}.bam', idx=read_files.get_indexes)
        get_final_bam_files
    output:
        bam = "{pfx}/{sample}/maps/final.bam"
    run:
        if len(input) > 1:
            shell('samtools merge -@4 {output.bam} {input}')
        else:
            # use hard link since input will be removed
            shell('ln {input} {output}')





# EOF
