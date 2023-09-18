
# utilities

rule index_csi:
    threads: 1
    input:
        "{pfx}.vcf.gz"
    output:
        "{pfx}.vcf.gz.csi"
    shell:
        "bcftools index --csi {input}"


rule index_tbi:
    threads: 1
    input:
        "{pfx}.vcf.gz"
    output:
        "{pfx}.vcf.gz.tbi"
    shell:
        "bcftools index --tbi {input}"

rule index_bai:
    threads: 1
    input:
        "{pfx}.bam"
    output:
        "{pfx}.bam.bai"
    shell:
        "samtools index {input}"
