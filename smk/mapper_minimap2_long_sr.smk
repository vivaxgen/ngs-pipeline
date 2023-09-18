
rule reads_mapping_long_sr:
    threads: 4
    input:
        read = "reads/{sample}.fastq.gz"
    output:
        bam = "maps/{sample}.bam"
    shell:
        "minimap2 -a {ref.mmi} -ax map-ont {input.read} | samtools sort -o {output.bam}"

# EOF
