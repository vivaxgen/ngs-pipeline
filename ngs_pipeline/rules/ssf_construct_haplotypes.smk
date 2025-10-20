
# This is a haplotype-constructing workflow to be used with singleton long read data
# This workflow utilizes var_call_ont.smk


include: "var_call_ont.smk"

rule haplotypes_report:
    localrule: True
    input:
        "gvcf/variants.g.vcf.gz",
        'logs/mapped-final.stats.txt',
        'logs/mapped-final.depth-base.tsv.gz',
        'logs/stats.tsv',
        'logs/depths.png'


rule ssf_vcf_noindels:
    threads: 1
    input:
        "vcf/variants.vcf.gz",
    output:
        "vcf/variants-noindels.vcf.gz"
    shell:
        "bcftools view -o {output} -V indels -i 'QUAL>10' {input}"

rule set_vcf_csq:
    threads: 1
    input:
        "vcf/variants-noindels.vcf.gz",
    output:
        "vcf/variants-noindels-csq.vcf.gz",
    shell:
        "bcftools csq -l -f {refseq} --gff {gff_file} -o {output} {input}"


rule ssf_construct_haplotyes:
    threads: 1
    input:
        vcf = "vcf/variants-noindels-csq.vcf.gz",
        bam = "maps/mapped-final.bam",
    output:
        vcf = "vcf/phased.vcf.gz",
        table = "haplotype-table.tsv"
    log:
        log1 = "logs/haplotype-report.json",
        log2 = "logs/construct-haplotypes.log"
    shell:
        "ngs-pl construct-haplotypes"
        "  --min-match-length 2000 --outvcf {output.vcf}"
        "  --min-qual 10"
        "  --cds-only"
        "  --min-depth-ratio 0.5"
        "  --variant-list-file {input.vcf}"
        "  --outlog {log.log1}"
        "  --outfile {output.table}"
        "  {input.bam}"
        "  2> {log.log2}"

# EOF
