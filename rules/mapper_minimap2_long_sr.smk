
rule reads_mapping_lr:
    threads: 8
    input:
        "{pfx}/{sample_id}/trimmed_reads/trimmed.fastq.gz"
    output:
        "{pfx}/{sample_id}/maps/sorted.bam"
    params:
        rg = lambda w: f"-R @RG\\\\tID:{w.sample_id}\\\\tSM:{w.sample_id}\\\\tLB:LIB-{w.sample_id}",
    shell: 
        "minimap2 -a {refmap} {params.rg} {input} "
        "|  samtools sort -@4 -o {output} "
        "&& sleep 2 && samtools index {output}"

# EOF
