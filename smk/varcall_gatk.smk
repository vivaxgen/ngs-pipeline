

rule gatk_haplotypecaller:
    input:
        "maps/mapped-dedup.bam"
    output:
        "gvcf/{sample}.g.vcf"
    shell:
        "gatk HaplotypeCaller --native-pair-hmm-threads 8 -R {refseq} -I {input} -ERC GVCF -O {output}"

