# ssf_varcall_clair3.smk - ngs-pipeline rules
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2025, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"


rule ssf_varcall_clair3:
    # sleep is needed to ensure symbolic links is fully created
    threads: thread_allocations.get('variant_calling', 4)
    input:
        bam = "maps/mapped-final.bam",
        idx = "maps/mapped-final.bam.bai",
    output:
        vcf = "vcf/clair3/variants.vcf.gz",
        vcf_tbi = "vcf/clair3/variants.vcf.gz.tbi",
    log:
        log1 = "logs/clair3.log",
        log2 = "logs/clair3.err",
    params:
        sample = lambda w: sample if "sample" not in w else w.sample,
        model = config.get('clair3_model', 'ont'),
        platform = config.get('clair3_platform', 'ont'),
        flags = config.get('clair3_flags', ''),
        extra_flags = config.get('clair3_extra_flags', ''),
        model_path = config.get('clair3_model_path', ''),
        tmpdir = config.get('tmpdir', '/tmp'),
        outdir = subpath(output.vcf, parent=True),
        outfmt = "",
    shell:
        "run_clair3.sh"
        "  --bam_fn {input.bam}"
        "  --ref_fn {refseq}"
        "  --threads 2 --platform {params.platform}"
        "  --output {params.outdir}"
        "  --model_path {params.model_path}"
        "  --no_phasing_for_fa"
        "  --sample_name={params.sample}"
        "  {params.outfmt}"
        "  {params.flags}"
        "  {params.extra_flags}"
        "  1> {log.log1} 2> {log.log2}"


rule ssf_varcall_symlink:
    localrule: True
    input:
        vcf = "vcf/clair3/merge_output.vcf.gz",
        tbi = "vcf/clair3/merge_output.vcf.gz.tbi",
    output:
        vcf = "vcf/variants.vcf.gz",
        tbi = "vcf/variants.vcf.gz.tbi",
    shell:
        "sleep 2"
        " && ln -srf {input.vcf} {output.vcf}"
        " && ln -srf {input.tbi} {output.tbi}"


use rule ssf_varcall_clair3 as ssf_varcall_clair3_gvcf with:
    output:
        vcf = "gvcf/clair3/merge_output.gvcf.gz",
        idx = "gvcf/clair3/merge_output.gvcf.gz.tbi",
    params:
        outfmt = "--gvcf",




rule ssf_varcall_fix_gvcf:
    input:
        gvcf = "gvcf/clair3/merge_output.gvcf.gz",
        tbi = "gvcf/clair3/merge_output.gvcf.gz.tbi",
    output:
        gvcf = f"gvcf/merge_output.fixed.g.vcf.gz",
        tbi = f"gvcf/merge_output.fixed.g.vcf.gz.tbi",
    shell:
        "gatk UpdateVcfSequenceDictionary"
        " -SD {refseq}"
        " -I {input.gvcf}"
        " -O {output.gvcf}"
        " && bcftools index -t {output.gvcf}"

rule ssf_varcall_sort_gvcf:
    localrule: True
    input:
        gvcf = "gvcf/merge_output.fixed.g.vcf.gz",
        tbi = "gvcf/merge_output.fixed.g.vcf.gz.tbi",
    output:
        gvcf = "gvcf/merge_output.sorted.g.vcf.gz",
        tbi = "gvcf/merge_output.sorted.g.vcf.gz.tbi",
    shell:
        "bcftools sort -o {output.gvcf} {input.gvcf}"
        " && bcftools index -t {output.gvcf}"


rule ssf_varcall_reformat_gvcf:
    input:
        gvcf = "gvcf/merge_output.sorted.g.vcf.gz",
        tbi = "gvcf/merge_output.sorted.g.vcf.gz.tbi",
    output:
        gvcf = f"gvcf/{sample}-{complete_region}.g.vcf.gz",
        tbi = f"gvcf/{sample}-{complete_region}.g.vcf.gz.tbi",
    shell:
        "bcftools annotate -x FORMAT/AF -o {output.gvcf} {input.gvcf}"
        " && bcftools index -t {output.gvcf}"


# EOF
