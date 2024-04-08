


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
        flags = config.get('baserecalibrator_flags', ''),
        extra_flags = config.get('baserecalibrator_extra_flags', ''),
    shell:
        "gatk {java_opts} BaseRecalibrator {params.flags} {params.extra_flags} "
        "-R {refseq} {params.known} {params.region_opts} -I {input.bam} -O {output} 2>{log}"


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
        flags = config.get('gatherbqsrr_flags', ''),
        extra_flags = config.get('gatherbqsrr_extra_flags', ''),
        n_input = lambda wildcards, input: len(input),
        input_files = lambda wildcards, input: '-I ' + ' -I '.join(input),
    shell:
        #input_opts = '-I ' + ' -I '.join(input)
        "gatk {java_opts} GatherBQSRReports {params.flags} {params.extra_flags} "
        "{params.input_files} -O {output} 2>{log}"


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
        flags = config.get('applybqsr_flags', ''),
        extra_flags = config.get('applybqsr_extra_flags', '')
    log:
        "logs/gatk-ApplyBQSR.log"
    shell:
        "gatk {java_opts} ApplyBQSR {params.flags} {params.extra_flags} "
        "-R {refseq} -I {input.bam} --bqsr-recal-file {input.table} -O {output} 2>{log}"

# EOF
