

rule gatk_haplotypecaller:
    threads: 2
    input:
        "maps/mapped-dedup-recal.bam"
    output:
        "gvcf/{sample}-{reg}.g.vcf.gz"
    shell:
        "gatk HaplotypeCaller --native-pair-hmm-threads 1 -R {refseq} -I {input} -L {wildcards.reg} -ploidy {ploidy} -ERC GVCF -O {output}"

