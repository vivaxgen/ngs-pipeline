__copyright__ = '''
run_joint_variant_caller.py - ngs-pipeline command line
[https://github.com/vivaxgen/ngs-pipeline]

(c) 2022-2023 Hidayat Trimarsanto <trimarsanto@gmail.com>

All right reserved.
This software is licensed under MIT license.
Please read the README.txt of this software.
'''

# to improve the responsiveness during bash autocomplete, do not import heavy
# modules (such as numpy, pandas, etc) here, but instead import them within the
# functions that requires respective heavy modules

import sys
import os
from ngs_pipeline import cerr, cexit, arg_parser, check_NGSENV_BASEDIR, snakeutils


def init_argparser():
    p = snakeutils.init_argparser(desc='run joint variant calling')
    p.arg_dict['snakefile'].choices = ['jointvarcall_gatk.smk', 'jointvarcall_freebayes.smk']

    # input/output options
    p.add_argument('-o', '--outdir', default='joint',
                   help='directory for output [joint/]')
    p.add_argument('source_dirs', nargs='+',
                   help='source directories containing sample directories')

    return p


def run_joint_variant_caller(args):

    check_NGSENV_BASEDIR()

    import datetime
    import pathlib
    import yaml
    import snakemake
    from snakemake.io import glob_wildcards

    # get snakefile to run
    if args.snakefile is None:
        if 'JOINTCALL_SMK' in os.environ:
            args.snakefile = os.environ['JOINTCALL_SMK']
        else:
            args.snakefile = 'jointvarcall_gatk.smk'
    cerr(f'Snakefile to be run: {args.snakefile}')

    # sanity check for all directory and duplicated sample codes
    source_dirs = [srcdir.removesuffix('/') for srcdir in args.source_dirs]
    sample_dirs = []
    sample_codes = {}
    duplicated_sample_codes = {}
    for srcdir in source_dirs:
        if not pathlib.Path(srcdir).is_dir():
            cexit(f'ERROR: directory {srcdir} does not exist!')
        S, = glob_wildcards(srcdir + '/{sample,[\\w-]+}')
        # filter for non-sample directories/files
        S = [s for s in S if s != 'config.yaml']
        for sample in S:
            if sample in sample_codes:
                if sample in duplicated_sample_codes:
                    duplicated_sample_codes[sample].append(srcdir)
                else:
                    duplicated_sample_codes[sample] = [sample_codes[sample], srcdir]
            else:
                sample_codes[sample] = srcdir
        sample_dirs += [f'{srcdir}/{s}' for s in S]

    if any(duplicated_sample_codes):
        cerr('ERROR:')
        cerr('Found duplicated sample code(s):')
        cerr(yaml.dump(duplicated_sample_codes))
        cexit('Pleae remove or rename the directory of those sample(s) to remove the duplication')

    config = dict(srcdirs=source_dirs, destdir=args.outdir.removesuffix('/'))
    status, elapsed_time = snakeutils.run_snakefile(args, config=config)

    if not status:
        cerr('[WARNING: joint varian calling step did not successfully complete]')
    cerr(f'[Finish joint variant calling (time: {elapsed_time})]')


def main(args):
    run_joint_variant_caller(args)


# EOF
