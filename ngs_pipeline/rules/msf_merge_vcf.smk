# msf_mapper_vcf.smk - ngs-pipeline rules
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2023, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"


rule merge_vcfs:
    input:
        vcfs = expand(f'{outdir}/{{sample}}/vcfs/variants.vcf.gz',
                      sample=read_files.samples()),
        idx = expand(f'{outdir}/{{sample}}/vcfs/variants.vcf.gz.csi',
                     sample=read_files.samples()),

    output:
        vcf = '{outdir}/merged.vcf.gz'
    shell:
        "bcftools merge -o {output.vcf} {input.vcfs}"


# EOF
