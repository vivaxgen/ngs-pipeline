# multistep_variant_caller.smk - ngs-pipeline rules
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2024 Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

outdir = config['outdir']
infiles = config['infiles']
underscore = config.get('underscore', 0)
singleton = config.get('singleton', False)
paired_end = config.get('paired_end', False)

include: 'utilities.smk'


rule all:
    input:
        f'{outdir}/concatenated.vcf.gz.tbi',
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
        infiles = ' '.join(infiles)
    shell:
        'ngs-pl generate-manifest -o {output} {params.underscore} --pause 3 '
        '{params.singleton} {params.paired_end} --pause 3 {infiles}'


rule run_prepare_sample_directory:
    localrule: True
    input:
        f'{outdir}/metafile/manifest.tsv'
    output:
        f'{outdir}/analysis/._prepared_'
    params:
        extra_opts = '--force'
    shell:
        'ngs-pl prepare-sample-directory {params.extra_opts} '
        '-o {outdir}/analysis -i {input} . '
        '&& touch {output}'


rule run_sample_variant_caller:
    localrule: True
    input:
        f'{outdir}/analysis/._prepared_'
    output:
        f'{outdir}/analysis/._completed_'
    shell:
        'ngs-pl run-sample-variant-caller --force --no-config-cascade '
        '--target all {outdir}/analysis '
        '&& touch {output}'

rule run_check_sample_variant_result:
    localrule: True
    input:
        f'{outdir}/analysis/._completed_'
    output:
        f'{outdir}/analysis/._checked_',
        f'{outdir}/failed_samples/._completed_'
    shell:
        'ngs-pl move-failed-samples -o {outdir}/failed_samples {outdir}/analysis/ '
        '&& touch {output[0]} && touch {output[1]}'


rule run_joint_variant_caller:
    localrule: True
    input:
        f'{outdir}/analysis/._checked_'
    output:
        f'{outdir}/joint/._completed_'
    shell:
        'ngs-pl run-joint-variant-caller --force --no-config-cascade '
        '-o {outdir}/joint {outdir}/analysis/ '
        '&& touch {output}'


rule concat_vcfs:
    threads: 1
    input:
        f'{outdir}/joint/._completed_'
    output:
        f'{outdir}/concatenated.vcf.gz'
    shell:
        'bcftools concat -o {output} {outdir}/joint/vcfs/*.vcf.gz'


rule consolidate_depth_base:
    localrule: True
    input:
        f'{outdir}/analysis/._checked_'
    output:
        f'{outdir}/reports/depth-base/._completed_'
    shell:
        '(for sample_dir in {outdir}/analysis/*; do'
        '     ln -sr ${{sample_dir}}/logs/mapped-final.depth-base.tsv.gz '
        '     {outdir}/reports/depth-base/`basename ${{sample_dir}}`.depth-base.tsv.gz;'
        ' done;) '
        '&& touch {output}'


rule consolidate_maps:
    localrule: True
    input:
        f'{outdir}/analysis/._checked_'
    output:
        f'{outdir}/reports/maps/._completed_'
    shell:
        '(for sample_dir in {outdir}/analysis/*; do'
        '     ln -sr ${{sample_dir}}/maps/mapped-final.bam '
        '     {outdir}/reports/maps/`basename ${{sample_dir}}`.bam;'
        ' done;) '
        '&& touch {output}'

# EOF
