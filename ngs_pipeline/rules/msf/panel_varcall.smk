# SPDX-FileCopyrightText: 2026 Hidayat Trimarsanto <trimarsanto@gmail.com>
# SPDX-License-Identifier: MIT

__copyright__ = "(c) 2026 Hidayat (Anto) Trimarsanto <trimarsanto@gmail.com>"
__author__ = "trimarsanto@gmail.com"
__license__ = "MIT"

from ngs_pipeline.rules import pkg

include: pkg("ngs_pipeline::msf/init.smk")
include: pkg(config.get("map_handler_wf", "ngs_pipeline::helper/map_handler.smk"))
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

if config.get("amplicon_based", False):
    amplicon_bed = get_abspath(config.get("amplicon_bed"))
    import pandas as pd
    
    amplicons_bed = pd.read_table(amplicon_bed, header=None)
    amplicons = amplicons_bed[3].unique().tolist()
    target_variants_vcf = dict()
    for amplicon in amplicons:
        if vcf_path := config.get(f"target_variants_vcf_{amplicon}", None):
            target_variants_vcf[amplicon] = get_abspath(vcf_path)


    print(target_variants_vcf.keys())
    wildcard_constraints:
        amplicon = "|".join(target_variants_vcf.keys())

    rule varcall_amplicon_concat:
        # note: the per-amplicon final vcf is deliberately NOT placed under a
        # "vcfs/" subdirectory (it lives at ".../{sample}/{amplicon}/variants.vcf.gz").
        # ".../{sample}/{amplicon}/vcfs/variants.vcf.gz" would collide with the
        # generic "{pfx}/{sample}/vcfs/variants.vcf.gz" pattern used by
        # rename_set_GT/clair3_symlink (pfx is an unconstrained wildcard), which
        # makes Snakemake treat {amplicon} as a bogus sample. Unlike the plain
        # "maps/mapped-final.bam" collision, this one isn't reliably resolved by
        # ruleorder: Snakemake's separate "no two rules make the same output"
        # consistency check (dag.py check_jobs) does not consult ruleorder at all
        # and can still raise AmbiguousRuleException depending on which path first
        # discovers the file.
        input:
            vcf = expand(f"{outdir}/samples/{{{{sample}}}}/{{amplicon}}/variants.vcf.gz", amplicon=target_variants_vcf.keys()),
            vcf_idx = expand(f"{outdir}/samples/{{{{sample}}}}/{{amplicon}}/variants.vcf.gz.csi", amplicon=target_variants_vcf.keys()),
        output:
            vcf = f"{outdir}/samples/{{sample}}/vcfs/variants.vcf.gz",
        shell:
            "bcftools concat -a {input.vcf} | bcftools sort -o {output.vcf}"
    
    rule split_according_to_amplicons:
        # note: the per-amplicon bam is deliberately NOT placed under a "maps/"
        # subdirectory. ".../{sample}/{amplicon}/maps/mapped-final.bam" would
        # collide with the generic "{pfx}/{sample}/maps/mapped-final.bam" pattern
        # used throughout map_handler.smk/clair3.smk (pfx is an unconstrained
        # wildcard), which makes Snakemake treat {amplicon} as a bogus sample and
        # silently drops this whole branch from the DAG when that lookup fails.
        input:
            bam_file = f"{outdir}/samples/{{sample}}/maps/mapped-final.bam",
            bai_file = f"{outdir}/samples/{{sample}}/maps/mapped-final.bam.bai",
        output:
            splitted_bam = expand(f"{outdir}/samples/{{{{sample}}}}/{{amplicon}}/mapped-final.bam", amplicon=target_variants_vcf.keys()),
            depths = expand(f"{outdir}/samples/{{{{sample}}}}/{{amplicon}}/depth.txt", amplicon=target_variants_vcf.keys()),
        params:
            output_dir = f"{outdir}/samples/{{sample}}",
            tolerate_start_bp = config.get("tolerate_start_bp", 5),
            tolerate_end_bp = config.get("tolerate_end_bp", 5),
        run:
            for _, row in amplicons_bed.iterrows():
                amplicon = row[3]
                if not amplicon in target_variants_vcf.keys():
                    continue
                chrom = row[0]
                start = row[1]
                end = row[2]
                tolerate_start = start - params.tolerate_start_bp
                tolerate_end = end + params.tolerate_end_bp
                outfile = f"{params.output_dir}/{amplicon}/mapped-final.bam"
                outdepth = f"{params.output_dir}/{amplicon}/depth.txt"
                shell("mkdir -p {params.output_dir}/{amplicon}")
                shell(f"samtools view -b {input.bam_file} -e 'pos >= {tolerate_start} && endpos <= {tolerate_end}' -o {outfile} {chrom}:{start}-{end}")
                shell(f"samtools index {outfile}")
                shell(f"samtools coverage -r {chrom}:{start}-{end} {outfile} > {outdepth}")

    match _varcaller:
        case "freebayes":
            use rule freebayes as freebayes_amplicon with:
                input:
                    bam = f"{outdir}/samples/{{sample}}/{{amplicon}}/mapped-final.bam",
                    idx = f"{outdir}/samples/{{sample}}/{{amplicon}}/mapped-final.bam.bai"
                output:
                    vcf = f"{outdir}/samples/{{sample}}/{{amplicon}}/variants.vcf.gz",
                params:
                    depth_file = f"{outdir}/samples/{{sample}}/{{amplicon}}/maps/depth.txt",
                    target = "",
                    vcf_target = lambda w: f"-@ {target_variants_vcf[w.amplicon]} "

            ruleorder: freebayes_amplicon > varcall_amplicon_concat > freebayes 

        case "clair3":
            use rule clair3 as clair3_amplicon with:
                input:
                    bam = f"{outdir}/samples/{{sample}}/{{amplicon}}/mapped-final.bam",
                    idx = f"{outdir}/samples/{{sample}}/{{amplicon}}/mapped-final.bam.bai",
                    model = f"{outdir}/samples/{{sample}}/reads/model-0.txt",
                output:
                    vcf = f"{outdir}/samples/{{sample}}/{{amplicon}}/clair3/merge_output.vcf.gz",
                    idx = f"{outdir}/samples/{{sample}}/{{amplicon}}/clair3/merge_output.vcf.gz.tbi",
                log:
                    log1 = f"{outdir}/samples/{{sample}}/{{amplicon}}/logs/clair3.log",
                    log2 = f"{outdir}/samples/{{sample}}/{{amplicon}}/logs/clair3.err",
                params:
                    vcf_target = lambda w: f' --vcf_fn={target_variants_vcf[w.amplicon]}' if target_variants_vcf else '',
            use rule rename_set_GT as rename_set_GT_amplicon with:
                input:
                    f"{outdir}/samples/{{sample}}/{{amplicon}}/clair3/merge_output.vcf.gz",
                output:
                    final = f"{outdir}/samples/{{sample}}/{{amplicon}}/variants.vcf.gz",

            ruleorder: rename_set_GT_amplicon > rename_set_GT > clair3_symlink
            ruleorder: clair3_amplicon > varcall_amplicon_concat > clair3
            ruleorder: varcall_amplicon_concat > rename_set_GT > clair3_symlink
        case _:
            raise ValueError(f"Unsupported varcaller_wf {config['varcaller_wf']} for amplicon-based panel variant calling")

    ruleorder: varcall_amplicon_concat > varcall

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