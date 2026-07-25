# SPDX-FileCopyrightText: 2023-2026 Hidayat Trimarsanto <trimarsanto@gmail.com>
# SPDX-License-Identifier: MIT

__copyright__ = "(c) 2023-2026 Hidayat (Anto) Trimarsanto <trimarsanto@gmail.com>"
__author__ = "trimarsanto@gmail.com"
__license__ = "MIT"

# set this up so joint variant caller knows which variant caller is used in this workflow
if "varcaller" in locals() and varcaller:
    cexit(f"varcaller is already defined as {varcaller}, cannot redefine it in freebayes.smk")
_varcaller = "freebayes"

# varcall_freebayes.smk - ngs-pipeline rules
# [https://github.com/vivaxgen/ngs-pipeline]

# targeted variant calling with freebayes, for either panel or discovery setting

# required variables:
# - refseq
# - target_variants (if panel variant calling)
# - min_read_qual

# optional config keys
# - vcf_variants
# - freebayes_extra_flags

vcf_variants = get_abspath(config["vcf_variants"]) if "vcf_variants" in config else ""

rule freebayes:
    threads: 2
    input:
        bam = "<sp>maps/mapped-final.bam",
        idx = "<sp>maps/mapped-final.bam.bai"
    output:
        vcf = "<sp>vcfs/variants.vcf.gz",
    params:
        target = f"--target {target_variants}" if target_variants else "",
        vcf_target = f"-@ {vcf_variants} -l" if vcf_variants else "",
        monomorphic = '--report-monomorphic' if target_variants else '',
        freebayes_extra_flags = config.get('freebayes_extra_flags', ''),
        min_read_qual = min_read_qual,
    shell:
        "freebayes -f {refseq} {params.target} {params.monomorphic} --haplotype-length 0 "
        "--min-base-quality {params.min_read_qual} {params.freebayes_extra_flags} {input.bam} "
        "| bcftools sort -o {output.vcf}"


# EOF
