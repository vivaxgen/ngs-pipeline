
# targeted variant calling with freebayes, for either panel or discovery setting

# required variables:
# - refseq
# - target_variants (if panel variant calling)
# - min_read_qual

# optional config keys
# - freebayes_extra_flags

rule freebayes:
    threads: 2
    input:
        bam = "{pfx}/maps/sorted.bam",
        idx = "{pfx}/maps/sorted.bam.bai"
    output:
        vcf = "{pfx}/vcfs/variants.vcf.gz",
    params:
        target = f'--target {target_variants}' if target_variants else '',
        monomorphic = '--report-monomorphic' if target_variants else '',
        freebayes_extra_flags = config.get('freebayes_extra_flags', ''),
    shell:
        "freebayes -f {refseq} {params.target} {params.monomorphic} --haplotype-length 0 "
        "--min-base-quality {min_read_qual} {params.freebayes_extra_flags} {input.bam} "
        "| bcftools sort -o {output.vcf}"

# EOF
