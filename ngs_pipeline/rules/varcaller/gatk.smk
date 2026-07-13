
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
        "<sp>maps/mapped-final-recal.bam"
    output:
        replace_sp("<sp>gvcf/{sample}-{reg}.g.vcf.gz"),
    log:
        "<sp>logs/haplotypecaller-{sample}-{reg}.log"
    params:
        sample = get_sample,
        reg = get_haplotypecaller_region,
        flags = config.get('haplotypecaller_flags', ''),
        extra_flags = config.get('haplotypecaller_extra_flags', ''),
    shell:
        "gatk {java_opts} HaplotypeCaller  --native-pair-hmm-threads {threads}"
        "  -R {refseq}  -I {input} {params.reg}  -ploidy {ploidy}  -ERC GVCF"
        "  {params.flags} {params.extra_flags}  -O {output}  2> {log}"

# EOF
