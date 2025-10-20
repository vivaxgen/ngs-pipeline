
rule fb_varcall:
    threads: 1
    input:
        'maps/mapped-final.bam'
    output:
        f'gvcf/{sample}-{{reg}}.g.vcf.gz'
    log:
        fb = 'logs/freebayes-{reg}.log'
    params:
        reg = get_interval,
        fb_args = "--haplotype-length -1 --min-coverage 10 " + config["freebayes_extra_flags"]
    shell:
        "freebayes --fasta-reference {refseq} --ploidy {ploidy} --min-alternate-count 2 {params.reg} "
        "{params.fb_args} {input} 2> {log.fb} "
        "| bcftools view -o {output} "
