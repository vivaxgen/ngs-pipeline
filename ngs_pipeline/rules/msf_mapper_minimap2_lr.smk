# msf_mapper_minimap2_lr.smk - ngs-pipeline rules
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2023, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

# required variables:
# - refmap

rule msf_mapping_lr:
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

