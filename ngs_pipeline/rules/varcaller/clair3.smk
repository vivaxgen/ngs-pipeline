# ssf_varcall_clair3.smk - ngs-pipeline rules
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2025, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

# set this up so joint variant caller knows which variant caller is used in this workflow
if "varcaller" in locals() and varcaller:
    cexit(f"varcaller is already defined as {varcaller}, cannot redefine it in clair3.smk")
_varcaller = "clair3"

ruleorder: clair3 > index_tbi
ruleorder: clair3_gvcf > index_tbi
ruleorder: clair3_symlink > index_tbi
ruleorder: clair3_fix_contigs_header > index_tbi
ruleorder: clair3_fix_gvcf > index_tbi
ruleorder: clair3_sort_gvcf > index_tbi
ruleorder: clair3_reformat_gvcf > index_tbi
ruleorder: clair3_fix_gvcf > clair3_reformat_gvcf


def get_model_flag(w, input):
    # read model name from the input.model

    if not input.model:
        return ""

    with open(input.model) as f:
        model_name = f.read().strip()

    model_path = f"{os.environ['VVG_BASEDIR']}/opt/clair3-models/{model_name}"
    if not os.path.exists(model_path):
        raise ValueError(f"Clair3 model path {model_path} does not exist.")

    return f"--model_path {model_path}"


rule clair3:
    threads: thread_allocations.get('variant_calling', 4)
    input:
        bam = "<sp>maps/mapped-final.bam",
        idx = "<sp>maps/mapped-final.bam.bai",
        # we assume the model is identical for all index in a sample, hence use model-0.txt
        model = "<sp>reads/model-0.txt" if ngs_platform.upper() in ['ONT'] else [],
    output:
        vcf = "<sp>vcfs/clair3/merge_output.vcf.gz",
        idx = "<sp>vcfs/clair3/merge_output.vcf.gz.tbi",
    log:
        log1 = "<sp>logs/clair3.log",
        log2 = "<sp>logs/clair3.err",
    params:
        sample = get_sample,
        platform = ngs_platform.lower(),
        flags = config.get('clair3_flags', ''),
        vcf_target = f' --vcf_fn={target_variants_vcf}' if target_variants_vcf else '',
        extra_flags = config.get('clair3_extra_flags', ''),
        model_flag = get_model_flag,
        outdir = subpath(output.vcf, parent=True),
        outfmt = "",
        # generate-null-gvcf params
        dict_file = f"{refseq.removesuffix('.fasta')}.dict",
        contig = lambda w: f"--contig {w.region}" if w.region != complete_region else "",
    shell:
        "run_clair3"
        "  --bam_fn {input.bam}"
        "  --ref_fn {refseq}"
        "  --threads {threads} --platform {params.platform}"
        "  --output {params.outdir}"
        "  {params.model_flag}"
        "  --sample_name={params.sample}"
        "  {params.vcf_target}"
        "  {params.outfmt}"
        "  {params.flags}"
        "  {params.extra_flags}"
        "  1> {log.log1} 2> {log.log2}"
        "  || true;"
        "  "
        "if [ ! -s {output.vcf} ]; then"
        "  ngs-pl generate-null-gvcf -o {output.vcf} -s {params.sample}"
        "  {params.contig} --dict {params.dict_file};"
        "  [ -f {output.idx} ] || tabix -p vcf {output.vcf};"
        "fi"


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


use rule clair3 as clair3_gvcf_region with:
    output:
        vcf = "<sp>gvcf/clair3-{region}/merge_output.gvcf.gz",
        idx = "<sp>gvcf/clair3-{region}/merge_output.gvcf.gz.tbi",
    log:
        log1 = "<sp>logs/clair3-{region}.log",
        log2 = "<sp>logs/clair3-{region}.err",
    params:
        outfmt = lambda w: f"--gvcf --include_all_ctgs" if (w.region == complete_region) else f"--gvcf --ctg_name {w.region}",


rule clair3_fix_contigs_header:
    input:
        gvcf = "<sp>gvcf/clair3-{region}/merge_output.gvcf.gz",
        idx = "<sp>gvcf/clair3-{region}/merge_output.gvcf.gz.tbi",
    output:
        gvcf = replace_sp("<sp>gvcf/clair3-{region}/merge_output.fixed.g.vcf.gz"),
        tbi = replace_sp("<sp>gvcf/clair3-{region}/merge_output.fixed.g.vcf.gz.tbi"),
    log:
        log1 = "<sp>logs/clair3_fix_contigs_header-{region}.log",
        log2 = "<sp>logs/clair3_fix_contigs_header-{region}.err",
    params:
        sample = get_sample,
        dict_file = f"{refseq.removesuffix('.fasta')}.dict",
    shell:
        "picard UpdateVcfSequenceDictionary"
        "  -I {input.gvcf}  -O {output.gvcf}"
        "  -SD {params.dict_file}"
        "  -CREATE_INDEX true"
        "  1> {log.log1} 2> {log.log2}"


rule clair3_fix_gvcf:
    input:
        gvcf = "<sp>gvcf/clair3-{region}/merge_output.fixed.g.vcf.gz",
        tbi = "<sp>gvcf/clair3-{region}/merge_output.fixed.g.vcf.gz.tbi",
    output:
        gvcf = replace_sp("<sp>gvcf/{sample}-{region}.g.vcf.gz"),
        tbi = replace_sp("<sp>gvcf/{sample}-{region}.g.vcf.gz.tbi"),
    log:
        log1 = "<sp>logs/clair3_fix_gvcf-{region}.log",
        log2 = "<sp>logs/clair3_fix_gvcf-{region}.err",
    params:
        sample = get_sample,
        dict_file = f"{refseq.removesuffix('.fasta')}.dict",
    shell:
        #"gatk UpdateVcfSequenceDictionary"
        #"bcftools reheader"
        #" --fai {refseq}.fai"
        #" -o {output.gvcf}"
        #" {input.gvcf}"
        #" && bcftools index -t {output.gvcf}"
        "picard SortVcf"
        "  -I {input.gvcf}  -O {output.gvcf}"
        "  -CREATE_INDEX true"
        "  1> {log.log1} 2> {log.log2}"


rule clair3_sort_gvcf:
    localrule: True
    input:
        gvcf = "{pfx}/merge_output.fixed.g.vcf.gz",
        tbi = "{pfx}/merge_output.fixed.g.vcf.gz.tbi",
    output:
        gvcf = "{pfx}/merge_output.sorted.g.vcf.gz",
        tbi = "{pfx}/merge_output.sorted.g.vcf.gz.tbi",
    shell:
        "bcftools sort -o {output.gvcf} {input.gvcf}"
        " && bcftools index -t {output.gvcf}"


rule clair3_reformat_gvcf:
    input:
        gvcf = "<sp>gvcf/clair3-{region}/merge_output.fixed.g.vcf.gz",
        tbi = "<sp>gvcf/clair3-{region}/merge_output.fixed.g.vcf.gz.tbi",
    output:
        gvcf = replace_sp("<sp>gvcf/{sample}-{region}.g.vcf.gz"),
        tbi = replace_sp("<sp>gvcf/{sample}-{region}.g.vcf.gz.tbi"),
    shell:
        "bcftools annotate -x FORMAT/AF -o {output.gvcf} {input.gvcf}"
        " && bcftools index -t {output.gvcf}"


# EOF
