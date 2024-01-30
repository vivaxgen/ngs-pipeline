
# required variables:
# - refmap

rule tgt_mapping_lr:
    threads: 8
    input:
        read = "{pfx}/{sample}/trimmed_reads/trimmed-{idx}.fastq.gz"
    output:
        bam = temp("{pfx}/{sample}/maps/sorted-{idx}.bam")
    params:
        rg = lambda w: f"-R @RG\\\\tID:{w.sample}-{w.idx}\\\\tSM:{w.sample}\\\\tLB:LIB-{w.sample}-{w.idx}",
        threads = lambda wildcards, threads: threads - 1,
    shell:
        "minimap2 -t {params.threads} -a {refmap} {params.rg} {input.read} "
        "| samtools sort -@4 -o {output.bam} "


# EOF

