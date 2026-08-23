
# set this up so joint variant caller knows which variant caller is used in this workflow
if "varcaller" in locals() and varcaller:
    cexit(f"varcaller is already defined as {varcaller}, cannot redefine it in gatk_drag.smk")
_varcaller = "gatk"

gatk_calibrate_str = config.get('gatk_calibrate_str', False)


def get_haplotypecaller_region(wildcards):
    if wildcards.reg == complete_region:
        if targetregion_file:
            return f'-L {targetregion_file}'
        return ''
    return f'-L {wildcards.reg}'


rule gatk_calibrate_STR:
    threads: 1
    input:
        bam = "<sp>maps/mapped-final.bam",
        bam_index = "<sp>maps/mapped-final.bam.bai",
    output:
        model = "<sp>maps/dragstr_model.txt"
    shell:
        "gatk {java_opts} CalibrateDragstrModel  -R {refseq}  -str {strtable_file}"
        "  -I {input.bam}"
        "  -O {output.model}"


rule gatk_drag_haplotypecaller:
    threads: thread_allocations.get('haplotyping', 2)
    input:
        # GATK DRAGEN use non-calibrated bam input
        #bam = "<sp>maps/mapped-final.bam",
        #bam_index = "<sp>maps/mapped-final.bam.bai",
        #model = "<sp>maps/dragstr_model.txt" if gatk_calibrate_str else [],
        bam = "{anypath}maps/mapped-final.bam",
        bam_index = "{anypath}maps/mapped-final.bam.bai",
        model = "{anypath}maps/dragstr_model.txt" if gatk_calibrate_str else [],
    output:
        #gvcf = replace_sp("<sp>gvcf/{sample}-{reg}.g.vcf.gz"),
        gvcf = "{anypath}gvcf/variants-{reg}.g.vcf.gz",
    log:
        "{anypath}logs/haplotypecaller-{reg}.log"
    params:
        sample = sample,
        reg = get_haplotypecaller_region,
        str_model = lambda w, input: f"--dragstr-params-path {input.model}" if gatk_calibrate_str else "",
        flags = config.get('haplotypecaller_flags', ''),
        extra_flags = config.get('haplotypecaller_extra_flags', ''),
    shell:
        "gatk {java_opts} HaplotypeCaller  --native-pair-hmm-threads {threads}"
        "  --dragen-mode true  {params.str_model}"
        "  -R {refseq}  -I {input.bam} {params.reg}  -ploidy {ploidy}  -ERC GVCF"
        "  {params.flags}  {params.extra_flags}  -O {output.gvcf} 2> {log}"


rule gatk_haplotypecaller_rename:
    # this is needed to make sure that the gvcf files have unique filenames
    # by having sample name, instead of just the region name
    localrule: True
    input:
        gvcf = "<sp>gvcf/variants-{reg}.g.vcf.gz",
        gvcf_index = "<sp>gvcf/variants-{reg}.g.vcf.gz.tbi"
    output:
        gvcf = replace_sp("<sp>gvcf/{sample}-{reg}.g.vcf.gz"),
        gvcf_index = replace_sp("<sp>gvcf/{sample}-{reg}.g.vcf.gz.tbi")
    shell:
        "ln {input.gvcf} {output.gvcf}"
        " && "
        "ln {input.gvcf_index} {output.gvcf_index}"


# EOF
