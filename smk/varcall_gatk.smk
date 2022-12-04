

rule gatk_haplotypecaller:
    threads: 1
    input:
        "maps/mapped-dedup-recal.bam"
    output:
        "gvcf/{sample}-{reg}.g.vcf"
    shell:
        "gatk HaplotypeCaller --native-pair-hmm-threads 1 -R {refseq} -I {input} -L {wildcards.reg} -ERC GVCF -O {output}"

