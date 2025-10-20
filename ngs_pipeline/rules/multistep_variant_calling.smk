# multistep_variant_caller.smk - ngs-pipeline rules
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2024 Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"


import pathlib

# -- preserved  keywords from snakeutils --

# config files to be passed down to any snakefile executor
configfiles = config['__configfiles__']
configfiles_args = ' '.join([f'-c {x}' for x in configfiles])


# -- input/output arguments --
outdir = config['outdir']
infiles = config['infiles']
underscore = config.get('underscore', 0)
singleton = config.get('singleton', False)
paired_end = config.get('paired_end', False)
manifest = config.get('manifest', None)
jobs = config.get('jobs', 32)
procfile = config.get('procfile', None)
rerun = config.get('rerun', False)
unlock = config.get('unlock', False)

# -- extra flags for each steps --
prepare_sample_directory_flags = config.get('prepare_sample_directory_flags', '')
sample_variant_caller_flags = config.get('sample_variant_caller_flags', '')
joint_variant_caller_flags = config.get('joint_variant_caller_flags', '')

# -- specific targets and rules --
sample_variant_caller_target = config.get('sample_variant_caller_target', 'all')
sample_variant_caller_wf = config.get("sample_variant_caller_wf", "var_call.smk")
joint_variant_caller_target = config.get('joint_variant_caller_target', 'all')
joint_variant_caller_wf = config.get('joint_variant_caller_wf', 'jointvarcall_gatk.smk')

include: "general_params.smk"
include: 'utilities.smk'


# retouch necessary touched files if rerun
if rerun or unlock:
    for fn in [f'{outdir}/analysis/._prepared_']:
        if (fp := pathlib.Path(fn)).exists():
            fp.touch()


rule all:
    input:
        f'{outdir}/joint/concatenated.vcf.gz.tbi',
        f'{outdir}/stats.tsv',
        f'{outdir}/reports/._completed_',


rule annotated_concatd_vcf:
    input:
        f'{outdir}/joint/concatenated-bcsq.vcf.gz.tbi',
        f'{outdir}/stats.tsv',
        f'{outdir}/reports/._completed_',


rule varcall_result:
    input:
        f'{outdir}/metafile/manifest.tsv',
        f'{outdir}/analysis/._completed_',
        f'{outdir}/joint/._completed_',


rule sample_variant_calling:
    input:
        branch(
            not unlock,
            then=[
                f'{outdir}/analysis/._completed_',
                f'{outdir}/completed_samples/._completed_',
                f'{outdir}/stats.tsv',
                f'{outdir}/reports/._completed_',
            ],
            otherwise=[
                f'{outdir}/analysis/._completed_',
            ]
        )


rule VCF:
    input:
        f'{outdir}/joint/._completed_',
        f'{outdir}/reports/._completed_',


rule concatenated_VCF:
    input:
        f'{outdir}/joint/concatenated.vcf.gz.tbi',
        f'{outdir}/reports/._completed_',


rule consolidated_reports:
    input:
        f'{outdir}/reports/depth-base/._completed_',
        f'{outdir}/reports/maps/._completed_',
    output:
        f'{outdir}/reports/._completed_'
    shell:
        'touch {output}'


rule generate_config:
    localrule: True
    output:
        config=f'{outdir}/metafile/config.yaml' if not rerun else [],
    params:
        config=config,
    run:
        import yaml

        if any(output.config):
            # remove any reserved keys (starting with __)
            config_copy = params.config.copy()
            for k in list(config_copy):
                if k.startswith('__'):
                    del config_copy[k]
            buff = yaml.safe_dump(config_copy)
            with open(output.config, 'w') as f:
                f.write("# this config yaml can still be overridden by command line\n" +
                        "# or by the config files nearest to working directory\n" +
                        "# if config_cascade is True\n\n")
                f.write(buff)


rule generate_manifest:
    localrule: True
    output:
        f'{outdir}/metafile/manifest.tsv' if not rerun else []
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
        flag = f'{outdir}/analysis/._prepared_',
        cfg = f'{outdir}/metafile/config.yaml',
    output:
        touch(f'{outdir}/analysis/._completed_')
    params:
        jobs = jobs,
        procfile = f'-P {procfile}' if procfile else '',
        rerun = '--rerun' if rerun else '',
        unlock = '--unlock' if unlock else '',
        target = sample_variant_caller_target,
        extra_flags = sample_variant_caller_flags,
    shell:
        'ngs-pl run-sample-variant-caller {params.extra_flags}'
        '  -c {input.cfg}'
        '  -j {params.jobs} {params.procfile} {params.rerun} {params.unlock}'
        '  --target {params.target}'
        '  --snakefile {sample_variant_caller_wf}'
        '  {outdir}/analysis'


rule run_check_sample_variant_result:
    localrule: True
    input:
        f'{outdir}/analysis/._completed_'
    output:
        touch(f'{outdir}/completed_samples/._completed_'),
        touch(f'{outdir}/failed_samples/._completed_'),
    shell:
        'ngs-pl consolidate-samples -o {outdir}/completed_samples '
        '-f {outdir}/failed_samples {outdir}/analysis/'


rule run_joint_variant_caller:
    localrule: True
    input:
        flag = f'{outdir}/completed_samples/._completed_',
        cfg = f'{outdir}/metafile/config.yaml',
    output:
        touch(f'{outdir}/joint/._completed_')
    params:
        rerun = '--rerun' if rerun else '',
        unlock = '--unlock' if unlock else '',
        extra_flags = joint_variant_caller_flags,
    shell:
        'ngs-pl run-joint-variant-caller {params.extra_flags}'
        '  {params.rerun} {params.unlock}'
        '  -c {input.cfg}'
        '  --snakefile {joint_variant_caller_wf}'
        '  --target {joint_variant_caller_target}'
        '  -o {outdir}/joint'
        '  {outdir}/completed_samples/'


rule concat_vcfs:
    threads: 1
    input:
        f'{outdir}/joint/._completed_'
    output:
        f'{outdir}/joint/concatenated.vcf.gz'
    shell:
        'bcftools concat -o {output} {outdir}/joint/vcfs/*.vcf.gz'

rule csq_vcf:
    input:
        vcf = "{fn}.vcf.gz",
        tbi = "{fn}.vcf.gz.tbi",
    output:
        vcf = "{fn}-bcsq.vcf.gz",
    shell:
        "bcftools csq -l -o {output.vcf} -f {refseq} --gff {gff_file} {input.vcf}"


rule gather_stats:
    localrule: True
    input:
        f'{outdir}/analysis/._completed_'
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
