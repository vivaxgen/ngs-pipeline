__copyright__ = '''
run_targeted_variant_caller.py - ngs-pipeline command line
[https://github.com/vivaxgen/ngs-pipeline]

(c) 2023 Hidayat Trimarsanto <trimarsanto@gmail.com>

All right reserved.
This software is licensed under MIT license.
Please read the README.txt of this software.
'''

# to improve the responsiveness during bash autocomplete, do not import heavy
# modules (such as numpy, pandas, etc) here, but instead import them within the
# functions that requires respective heavy modules

import os
from ngs_pipeline import (cerr, cexit, arg_parser,
                          check_NGSENV_BASEDIR, check_NGS_PIPELINE_BASE,
                          get_snakefile_path, setup_config, snakeutils)


def init_argparser():
    p = snakeutils.init_argparser(desc='run targeted variant calling')
    p.arg_dict['snakefile'].choices = [
        'panel_varcall_pe.smk',
        'panel_varcall_lr.smk',
        'msf_panel_varcall_pe.smk',
        'msf_panel_varcall_lr.smk'
    ]

    # input/output options
    p.add_argument('-u', '--underscore', default=0, type=int,
                   help='number of undercore character to be stripped, counted in reverse')

    p.add_argument('-o', '--outdir', default='analysis',
                   help='directory for output [analysis/]')
    p.add_argument('infiles', nargs='+',
                   help='FASTQ input files, eg. sample-1.fastq.gz')

    return p


def run_targeted_variant_caller(args):

    config = dict(
        infiles=args.infiles,
        underscore=args.underscore,
        outdir=args.outdir
    )
    status, elapsed_time = snakeutils.run_snakefile(args, config=config)

    if not status:
        cerr('[WARNING: targeted variant calling did not successfully complete]')
    cerr(f'[Finish targeted variant calling (time: {elapsed_time})]')


def main(args):
    run_targeted_variant_caller(args)

# EOF
