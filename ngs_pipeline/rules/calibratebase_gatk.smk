


rule gatk_baserecalibrator:
    # this rule generates calibration table from deduplicated bam file for each region
    threads: 1
    input:
        bam = 'maps/mapped-final.bam',
        # the following is for sanity check only
        known = f'{knownvariants_dir}/{{reg}}.bed.gz',
    output:
        temp("maps/recal-{reg}.table")
    log:
        "logs/gatk-BaseRecalibrator-{reg}.log"
    params:
        sample = sample,
        known = f"--known-sites {knownvariants_dir}/{{reg}}.bed.gz",
        region_opts = '-L {reg}',
    shell:
        "gatk BaseRecalibrator -R {refseq} {params.known} {params.region_opts} -I {input.bam} -O {output} 2>{log}"


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
        sample = sample,
        n_input = lambda wildcards, input: len(input),
    run:
        input_opts = '-I ' + ' -I '.join(input)
        shell("gatk {java_opts} GatherBQSRReports {input_opts} -O {output} 2>{log}")


rule gatk_applybqsr:
    # this rule applies calibration from the single calibration table
    threads: 2
    input:
        bam = "maps/mapped-final.bam",
        table = "maps/recal.table"
    output:
        "maps/mapped-final-recal.bam" if keep_recalibrated_bam else temp("maps/mapped-final-recal.bam")
    params:
        sample = sample,
    log:
        "logs/gatk-ApplyBQSR.log"
    shell:
        "gatk {java_opts} ApplyBQSR -R {refseq} -I {input.bam} --bqsr-recal-file {input.table} -O {output} 2>{log}"

# EOF
