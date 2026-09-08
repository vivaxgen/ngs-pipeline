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
# - target_variants_vcf
# - freebayes_extra_flags

target_variants_vcf = get_abspath(config["target_variants_vcf"]) if "target_variants_vcf" in config else ""

# https://github.com/freebayes/freebayes/issues/764 targeted varcall to handle low quality samples with low reads- freebayes does not fail gracefully with -@ -l
rule freebayes:
    threads: 2
    input:
        bam = "<sp>maps/mapped-final.bam",
        idx = "<sp>maps/mapped-final.bam.bai"
    output:
        vcf = "<sp>vcfs/variants.vcf.gz",
    params:
        sample = get_sample,
        depth_file = None,
        mindepth = config.get("min_depth", 10),
        target = f"--target {target_variants}" if target_variants else "",
        vcf_target = f"-@ {target_variants_vcf} " if target_variants_vcf else "",
        input_allele_only = "-l" if config.get("input_allele_only", False) else "",
        monomorphic = '--report-monomorphic' if target_variants else '',
        freebayes_extra_flags = config.get('freebayes_extra_flags', ''),
        min_read_qual = min_read_qual,
    run:
        import os
        import pandas as pd
        if params.depth_file:
            if os.path.exists(params.depth_file):
                df = pd.read_table(params.depth_file)
                depth = df["meandepth"].values[0]
                if depth < params.mindepth:
                    with open(output.vcf.replace(".vcf.gz", ".vcf"), "w") as f:
                        f.write(f"##fileformat=VCFv4.2\n")
                        f.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{params.sample}\n")
                    shell(f"bgzip {output.vcf.replace('.vcf.gz', '.vcf')}")
        if not os.path.exists(output.vcf):
            shell(
                "freebayes -f {refseq} {params.target} {params.vcf_target} {params.input_allele_only} {params.monomorphic} --haplotype-length 0 "
                "--min-base-quality {params.min_read_qual} {params.freebayes_extra_flags} {input.bam} "
                "| bcftools sort -o {output.vcf}"
            )

# EOF
