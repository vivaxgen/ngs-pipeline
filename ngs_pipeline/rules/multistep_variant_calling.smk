# multistep_variant_caller.smk - ngs-pipeline rules
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2024 Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

outdir = config['outdir']
infiles = config['infiles']
underscore = config.get('underscore', 0)
singleton = config.get('singleton', False)
paired_end = config.get('paired_end', False)
manifest = config.get('manifest', None)
jobs = config.get('jobs', 32)

prepare_sample_directory_flags = config.get('prepare_sample_directory_flags', '')
sample_variant_caller_flags = config.get('sample_variant_caller_flags', '')

include: 'utilities.smk'


rule all:
    input:
        f'{outdir}/concatenated.vcf.gz.tbi',
        f'{outdir}/stats.tsv',
        f'{outdir}/reports/._completed_',


rule varcall_result:
    input:
        f'{outdir}/metafile/manifest.tsv',
        f'{outdir}/analysis/._completed_',
        f'{outdir}/joint/._completed_'


rule concat_result:
    input:
        f'{outdir}/concatenated.vcf.gz.tbi'


rule consolidated_reports:
    input:
        f'{outdir}/reports/depth-base/._completed_',
        f'{outdir}/reports/maps/._completed_',
    output:
        f'{outdir}/reports/._completed_'
    shell:
        'touch {output}'


rule generate_manifest:
    localrule: True
    output:
        f'{outdir}/metafile/manifest.tsv'
    params:
        underscore = f'--underscore {underscore}' if underscore else '',
        singleton = '--single' if singleton else '',
        paired_end = '--paired' if paired_end else '',
        manifest = f'-i {manifest}' if manifest else '',
        infiles = ' '.join(infiles)
    shell:
        'ngs-pl generate-manifest -o {output} {params.underscore} --pause 5 '
        '{params.singleton} {params.paired_end} {params.manifest} {infiles}'


rule run_prepare_sample_directory:
    localrule: True
    input:
        f'{outdir}/metafile/manifest.tsv',
    output:
        touch(f'{outdir}/analysis/._prepared_')
    params:
        extra_flags = prepare_sample_directory_flags
    shell:
        'ngs-pl prepare-sample-directory {params.extra_flags} '
        '-o {outdir}/analysis -i {input} .'


rule run_sample_variant_caller:
    localrule: True
    input:
        f'{outdir}/analysis/._prepared_'
    output:
        touch(f'{outdir}/analysis/._completed_')
    params:
        jobs = jobs,
        extra_flags = sample_variant_caller_flags
    shell:
        'ngs-pl run-sample-variant-caller {params.extra_flags} '
        '-j {params.jobs} '
        '--target all {outdir}/analysis'


rule run_check_sample_variant_result:
    localrule: True
    input:
        f'{outdir}/analysis/._completed_'
    output:
        touch(f'{outdir}/completed_samples/._completed_'),
        touch(f'{outdir}/failed_samples/._completed_')
    shell:
        'ngs-pl consolidate-samples -o {outdir}/completed_samples '
        '-f {outdir}/failed_samples {outdir}/analysis/'


rule run_joint_variant_caller:
    localrule: True
    input:
        f'{outdir}/completed_samples/._completed_'
    output:
        touch(f'{outdir}/joint/._completed_')
    shell:
        'ngs-pl run-joint-variant-caller --force --no-config-cascade '
        '-o {outdir}/joint {outdir}/analysis/'


rule concat_vcfs:
    threads: 1
    input:
        f'{outdir}/joint/._completed_'
    output:
        f'{outdir}/concatenated.vcf.gz'
    shell:
        'bcftools concat -o {output} {outdir}/joint/vcfs/*.vcf.gz'


rule gather_stats:
    localrule: True
    input:
        f'{outdir}/completed_samples/._completed_'
    output:
        f'{outdir}/stats.tsv'
    shell:
        'ngs-pl gather-stats -o {output} {outdir}/analysis'


rule consolidate_depth_base:
    localrule: True
    input:
        f'{outdir}/completed_samples/._completed_'
    output:
        touch(f'{outdir}/reports/depth-base/._completed_')
    shell:
        'ngs-pl consolidate-reports -o {outdir}/reports/depth-base '
        '-t depth-base {outdir}/analysis/'


rule consolidate_maps:
    localrule: True
    input:
        f'{outdir}/completed_samples/._completed_'
    output:
        touch(f'{outdir}/reports/maps/._completed_')
    shell:
        'ngs-pl consolidate-reports -o {outdir}/reports/maps '
        '-t map {outdir}/analysis/'


# EOF
