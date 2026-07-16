# ssf_varcall_clair3.smk - ngs-pipeline rules
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2025, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"


ruleorder: clair3 > index_tbi

rule clair3:
    # sleep is needed to ensure symbolic links is fully created
    threads: thread_allocations.get('variant_calling', 4)
    input:
        bam = "<sp>maps/mapped-final.bam",
        idx = "<sp>maps/mapped-final.bam.bai",
    output:
        vcf = "<sp>vcfs/clair3/merge_output.vcf.gz",
        vcf_tbi = "<sp>vcfs/clair3/merge_output.vcf.gz.tbi",
    log:
        log1 = "<sp>logs/clair3.log",
        log2 = "<sp>logs/clair3.err",
    params:
        sample = get_sample,
        platform = config.get('clair3_platform', 'ont'),
        flags = config.get('clair3_flags', ''),
        vcf_target = f' --vcf_fn={target_variants_vcf}' if target_variants_vcf else '',
        extra_flags = config.get('clair3_extra_flags', ''),
        model_path = config.get('clair3_model_path', ''),
        outdir = subpath(output.vcf, parent=True),
        outfmt = "",
    shell:
        "run_clair3"
        "  --bam_fn {input.bam}"
        "  --ref_fn {refseq}"
        "  --threads {threads} --platform {params.platform}"
        "  --output {params.outdir}"
        "  --model_path {params.model_path}"
        "  --sample_name={params.sample}"
        "  {params.vcf_target}"
        "  {params.outfmt}"
        "  {params.flags}"
        "  {params.extra_flags}"
        "  1> {log.log1} 2> {log.log2}"


rule clair3_symlink:
    localrule: True
    input:
        vcf = "<sp>vcfs/clair3/merge_output.vcf.gz",
        tbi = "<sp>vcfs/clair3/merge_output.vcf.gz.tbi",
    output:
        vcf = "<sp>vcfs/variants.vcf.gz",
        tbi = "<sp>vcfs/variants.vcf.gz.tbi",
    shell:
        "sleep 2"
        " && ln -srf {input.vcf} {output.vcf}"
        " && ln -srf {input.tbi} {output.tbi}"


use rule clair3 as clair3_gvcf with:
    output:
        vcf = "<sp>gvcf/clair3/merge_output.gvcf.gz",
        idx = "<sp>gvcf/clair3/merge_output.gvcf.gz.tbi",
    params:
        outfmt = "--gvcf",


rule clair3_fix_gvcf:
    input:
        gvcf = "<sp>gvcf/clair3/merge_output.gvcf.gz",
        tbi = "<sp>gvcf/clair3/merge_output.gvcf.gz.tbi",
    output:
        gvcf = "<sp>gvcf/merge_output.fixed.g.vcf.gz",
        tbi = "<sp>gvcf/merge_output.fixed.g.vcf.gz.tbi",
    shell:
        "gatk UpdateVcfSequenceDictionary"
        " -SD {refseq}"
        " -I {input.gvcf}"
        " -O {output.gvcf}"
        " && bcftools index -t {output.gvcf}"

rule clair3_sort_gvcf:
    localrule: True
    input:
        gvcf = "<sp>gvcf/merge_output.fixed.g.vcf.gz",
        tbi = "<sp>gvcf/merge_output.fixed.g.vcf.gz.tbi",
    output:
        gvcf = "<sp>gvcf/merge_output.sorted.g.vcf.gz",
        tbi = "<sp>gvcf/merge_output.sorted.g.vcf.gz.tbi",
    shell:
        "bcftools sort -o {output.gvcf} {input.gvcf}"
        " && bcftools index -t {output.gvcf}"


rule clair3_reformat_gvcf:
    input:
        gvcf = "<sp>gvcf/merge_output.sorted.g.vcf.gz",
        tbi = "<sp>gvcf/merge_output.sorted.g.vcf.gz.tbi",
    output:
        gvcf = f"gvcf/{sample}-{complete_region}.g.vcf.gz",
        tbi = f"gvcf/{sample}-{complete_region}.g.vcf.gz.tbi",
    shell:
        "bcftools annotate -x FORMAT/AF -o {output.gvcf} {input.gvcf}"
        " && bcftools index -t {output.gvcf}"


# EOF
