
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
        bam = "maps/mapped-final.bam",
        bam_index = "maps/mapped-final.bam.bai",
    output:
        model = "maps/dragstr_model.txt"
    shell:
        "gatk {java_opts} CalibrateDragstrModel  -R {refseq}  -str {strtable_file}"
        "  -I {input.bam}"
        "  -O {output.model}"


rule gatk_drag_haplotypecaller:
    threads: thread_allocations.get('haplotyping', 2)
    input:
        # GATK DRAGEN use non-calibrated bam input
        bam = "maps/mapped-final.bam",
        bam_index = "maps/mapped-final.bam.bai",
        model = "maps/dragstr_model.txt" if gatk_calibrate_str else [],
    output:
        gvcf = "gvcf/{sample}-{reg}.g.vcf.gz",
    log:
        "logs/haplotypecaller-{sample}-{reg}.log"
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


# EOF
