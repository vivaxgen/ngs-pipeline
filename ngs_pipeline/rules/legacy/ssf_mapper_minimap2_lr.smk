# msf_mapper_minimap2_lr.smk - ngs-pipeline rules
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2025, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

# required variables:
# - refmap

rule ssf_mapping_lr:
    threads: 8
    input:
        read = "trimmed-reads/trimmed-{idx}_R0.fastq.gz"
    output:
        bam = temp("maps/mapped-sorted-{idx}.bam")
    log:
        log1 = "logs/minimap2-{idx}.log",
    params:
        sample = sample,
        rg = lambda w: f"-R @RG\\\\tID:{sample}-{w.idx}\\\\tSM:{sample}\\\\tLB:LIB-{sample}-{w.idx}",
        threads = lambda wildcards, threads: threads - 1,
        flags = config.get('minimap2_flags', ''),
        extra_flags = config.get('minimap2_extra_flags', ''),
    shell:
        "minimap2 -t {params.threads} --cs --ds --MD --eqx -a {refmap}"
        "  {params.flags} {params.rg} {params.extra_flags}"
        "  {input.read} 2> {log.log1}"
        " | samtools sort -@4 -o {output.bam}"

# EOF