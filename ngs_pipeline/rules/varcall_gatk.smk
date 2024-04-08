

rule gatk_haplotypecaller:
    threads: thread_allocations.get('haplotyping', 2)
    input:
        "maps/mapped-final-recal.bam"
    output:
        "gvcf/{sample}-{reg}.g.vcf.gz",
    log:
        "logs/haplotypecaller-{sample}-{reg}.log"
    params:
        sample = sample,
        flags = config.get('haplotypecaller_flags', ''),
        extra_flags = config.get('haplotypecaller_extra_flags', ''),
    shell:
        "gatk {java_opts} HaplotypeCaller --native-pair-hmm-threads 1 "
        "-R {refseq} -I {input} -L {wildcards.reg} -ploidy {ploidy} -ERC GVCF "
        "{params.flags} {params.extra_flags} -O {output} 2> {log}"

# EOF
