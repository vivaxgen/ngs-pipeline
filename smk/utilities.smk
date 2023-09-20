
# utilities

rule index_csi:
    threads: 1
    input:
        "{pfx}/{reg}.vcf.gz"
    output:
        "{pfx}/{reg}.vcf.gz.csi"
    shell:
        "bcftools index --csi {input}"


rule index_tbi:
    threads: 1
    input:
        "{pfx}/{reg}.vcf.gz"
    output:
        "{pfx}/{reg}.vcf.gz.tbi"
    shell:
        "bcftools index --tbi {input}"

rule index_bai:
    threads: 1
    input:
        "{pfx}/{reg}.bam"
    output:
        "{pfx}/{reg}.bam.bai"
    shell:
        "samtools index {input}"
