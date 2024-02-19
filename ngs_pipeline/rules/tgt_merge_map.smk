
# required params
# - read_files


def get_sorted_bam_files(w):
    return expand('{{pfx}}/{{sample}}/maps/sorted-{idx}.bam', idx=read_files.get_indexes(w.sample))

rule merge_map:
    # this rule merges dedup input BAM
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


# EOF