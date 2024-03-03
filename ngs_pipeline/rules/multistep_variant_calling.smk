# multistep_variant_caller.smk - ngs-pipeline rules
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2024 Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

outdir = config['outdir']
infiles = config['infiles']


rule all:
    input:
        f'{outdir}/metafile/manifest.tsv',
        f'{outdir}/analysis/._completed_',
        f'{outdir}/joint/._completed_'


rule generate_manifest:
    localrule: True
    output:
        f'{outdir}/metafile/manifest.tsv'
    params:
        underscore = config['underscore'],
        infiles = ' '.join(infiles)
    shell:
        'ngs-pl generate-manifest -o {output} -u {params.underscore} --pause 3 '
        '{infiles} '


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
        '-o {outdir}/analysis -i {input} .'
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

# EOF
