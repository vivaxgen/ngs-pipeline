# SPDX-FileCopyrightText: 2026 Hidayat Trimarsanto <trimarsanto@gmail.com>
# SPDX-License-Identifier: MIT

__copyright__ = "(c) 2026 Hidayat (Anto) Trimarsanto <trimarsanto@gmail.com>"
__author__ = "trimarsanto@gmail.com"
__license__ = "MIT"

from ngs_pipeline.rules import pkg

include: pkg("ngs_pipeline::msf/init.smk")
include: pkg("ngs_pipeline::helper/map_handler.smk")
include: pkg("ngs_pipeline::helper/genereport.smk")

# include trimmer
include: pkg(config["trimmer_wf"])

# include mapper
include: pkg(config["mapper_wf"])

# include variant caller
include: pkg(config["varcaller_wf"])

rule all:
    input:
        f"{outdir}/merged.vcf.gz",
        f"{outdir}/stats.tsv",


rule mapping:
    input:
        expand(f"{outdir}/samples/{{sample}}/maps/mapped-final.bam", sample=read_files.samples()),
        expand(f"{outdir}/samples/{{sample}}/logs/stats.tsv", sample=read_files.samples()),
        f"{outdir}/stats.tsv",


rule mapping_stats:
    localrule: True
    input:
        expand(f"{outdir}/samples/{{sample}}/logs/stats.tsv", sample=read_files.samples()),
    output:
        f"{outdir}/stats.tsv"
    shell:
        'ngs-pl gather-stats -o {output} {outdir}/samples'


rule varcall:
    input:
        expand(f"{outdir}/samples/{{sample}}/vcfs/variants.vcf.gz", sample=read_files.samples()),


rule merge_vcfs:
    input:
        vcfs = expand(f"{outdir}/samples/{{sample}}/vcfs/variants.vcf.gz",
                      sample=read_files.samples()),
        idx = expand(f"{outdir}/samples/{{sample}}/vcfs/variants.vcf.gz.csi",
                     sample=read_files.samples()),

    output:
        vcf = f"{outdir}/merged.vcf.gz"
    shell:
        "bcftools merge -m all -o {output.vcf} {input.vcfs}"


# EOF