# msf_mapper_minimap2_pe.smk - ngs-pipeline rules
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2023, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

# required variables:
# - refmap

rule msf_mapping:
    threads: 8
    input:
        read1 = "{pfx}/{sample}/trimmed_reads/trimmed-{idx}_R1.fastq.gz"
        read2 = "{pfx}/{sample}/trimmed_reads/trimmed-{idx}_R2.fastq.gz"
    output:
        bam = "{pfx}/{sample}/maps/sorted-{idx}.bam"
    params:
        rg = lambda w: f"-R @RG\\\\tID:{w.sample}-{w.idx}\\\\tSM:{w.sample}\\\\tLB:LIB-{w.sample}-{w.idx}",
        threads = lambda wildcards, threads: threads - 1,
    shell:
        "minimap2 -t {params.threads} -a {refmap} {params.rg} {input.read1} {input.read2} "
        "| samtools sort -@4 -o {output.bam} "

# EOF
