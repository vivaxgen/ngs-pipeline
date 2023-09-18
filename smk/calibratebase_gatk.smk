

rule gatk_baserecalibrator:
    # this rule generates calibration table from deduplicated bam file for each region
    threads: 1
    input:
        "maps/mapped-dedup.bam"
    output:
        temp("maps/recal-{reg}.table")
    params:
        known = f"--known-sites {knownsites_file}",
        region_opts = '-L {reg}'
    shell:
        "gatk BaseRecalibrator -R {refseq} {params.known} {params.region_opts} -I {input} -O {output}"


rule gatk_gatherbsqr:
    # this rule merges all calibration tables into single table
    threads: 2
    input:
        expand('maps/recal-{reg}.table', reg=REGIONS)
    output:
        "maps/recal.table"
    log:
        "logs/gatk-GatherBSQRReports.log"
    params:
        n_input = lambda wildcards, input: len(input)
    run:
        input_opts = '-I ' + ' -I '.join(input)
        shell("gatk {java_opts} GatherBQSRReports {input_opts} -O {output} 2>{log}")


rule gatk_applybqsr:
    # this rule applies calibration from the single calibration table
    threads: 2
    input:
        bam = "maps/mapped-dedup.bam",
        table = "maps/recal.table"
    output:
        "maps/mapped-dedup-recal.bam" if keep_recalibrated_bam else temp("maps/mapped-dedup-recal.bam")
    shell:
        "gatk {java_opts} ApplyBQSR -R {refseq} -I {input.bam} --bqsr-recal-file {input.table} -O {output}"

# EOF

