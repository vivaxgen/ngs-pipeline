
rule gatk_baserecalibrator_old:
    input:
        "maps/mapped-dedup.bam"
    output:
        "maps/mapped-dedup-recal.bam-XXX"
    params:
        known = f"--known-sites {knownsites_file}",
        table = "maps/recal.table"
    shell:
        """
        gatk BaseRecalibrator -R {refseq} {params.known} -I {input} -O {params.table}
        gatk {java_opts} ApplyBQSR -R {refseq} -I {input} --bqsr-recal-file {params.table} -O {output}
        """

rule gatk_baserecalibrator:
    threads: 1
    input:
        "maps/mapped-dedup.bam"
    output:
        "maps/recal-{reg}.table"
    params:
        known = f"--known-sites {knownsites_file}",
        region_opts = '-L {reg}'
    shell:
        "gatk BaseRecalibrator -R {refseq} {params.known} {params.region_opts} -I {input} -O {output}"


rule gatk_gatherbsqr:
    threads: 2
    input:
        expand('maps/recal-{reg}.table', reg=REGIONS)
    output:
        "maps/recal.table"
    run:
        input_opts = '-I ' + ' -I '.join(input)
        shell("gatk {java_opts} GatherBQSRReports {input_opts} -O {output}")


rule gatk_applybqsr:
    threads: 2
    input:
        bam = "maps/mapped-dedup.bam",
        table = "maps/recal.table"
    output:
        "maps/mapped-dedup-recal.bam"
    shell:
        "gatk {java_opts} ApplyBQSR -R {refseq} -I {input.bam} --bqsr-recal-file {input.table} -O {output}"

# EOF

