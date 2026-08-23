
# set this up so joint variant caller knows which variant caller is used in this workflow
if "varcaller" in locals() and varcaller:
    cexit(f"varcaller is already defined as {varcaller}, cannot redefine it in gatk.smk")
_varcaller = "gatk"

include: config.get("base_calibrator_wf", "gatk_calibratebase.smk")

def get_haplotypecaller_region(wildcards):
    if wildcards.reg == complete_region:
        if targetregion_file:
            return f'-L {targetregion_file}'
        return ''
    return f'-L {wildcards.reg}'


rule gatk_haplotypecaller:
    threads: thread_allocations.get('haplotyping', 2)
    input:
        #"<sp>maps/mapped-final-recal.bam"
        "{anypath}maps/mapped-final.bam"
    output:
        "{anypath}gvcf/variants-{reg}.g.vcf.gz",
    log:
        "{anypath}logs/haplotypecaller-{reg}.log"
    params:
        sample = get_sample,
        reg = get_haplotypecaller_region,
        flags = config.get('haplotypecaller_flags', ''),
        extra_flags = config.get('haplotypecaller_extra_flags', ''),
    shell:
        "gatk {java_opts} HaplotypeCaller  --native-pair-hmm-threads {threads}"
        "  -R {refseq}  -I {input} {params.reg}  -ploidy {ploidy}  -ERC GVCF"
        "  {params.flags} {params.extra_flags}  -O {output}  2> {log}"


rule gatk_haplotypecaller_rename:
    # this is needed to make sure that the gvcf files have unique filenames
    # by having sample name, instead of just the region name
    localrule: True
    input:
        gvcf = "<sp>gvcf/variants-{reg}.g.vcf.gz",
        gvcf_index = "<sp>gvcf/variants-{reg}.g.vcf.gz.tbi",
    output:
        gvcf = replace_sp("<sp>gvcf/{sample}-{reg}.g.vcf.gz"),
        gvcf_index = replace_sp("<sp>gvcf/{sample}-{reg}.g.vcf.gz.tbi"),
    shell:
        "ln {input.gvcf} {output.gvcf}"
        " && "
        "ln {input.gvcf_index} {output.gvcf_index}"


# EOF
