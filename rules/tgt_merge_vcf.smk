

rule merge_vcfs:
    input:
        expand()
    output:
        vcf = '{outdir}/merged.vcf.gz'
    shell:
        "bcftools merge -o {output.vcf} {input}"


# EOF
