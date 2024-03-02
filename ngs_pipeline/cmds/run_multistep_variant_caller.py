# run_multistep_variant_caller.py - ngs-pipeline command
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

# to improve the responsiveness during bash autocomplete, do not import heavy
# modules (such as numpy, pandas, etc) here, but instead import them within the
# functions that require the respective heavy modules

import os
import pathlib
from ngs_pipeline import (cerr, cexit, arg_parser, get_snakefile_path,
                          snakeutils)

# this is a wrapper to run sample_preparation, individual sample genotyping,
# and joint variant calling in a single command


def init_argparser():
    p = snakeutils.init_argparser('run whole steps of discovery variant caller')
    p.add_argument('-u', '--underscore', type=int, default=0,
                   help='no of consecutive underscore to be stripped from '
                   'filenames to form sample code, counted in reverse')
    
    p.add_argument('-o', '--outdir',
                   help='output directory')
    p.add_argument('infiles', nargs='+',
                   help="read files")
    
    return p


def run_multistep_variant_valler(args):

    import ngs_pipeline
    from ngs_pipeline import get_snakefile_path, snakeutils
    
    config=dict(
        underscore=args.underscore,
        outdir=args.outdir,
        infiles=args.infiles,
    )

    args.snakefile = get_snakefile_path(
        'multistep_variant_calling.smk',
        from_module=ngs_pipeline
    )

    status, elapsed_time = snakeutils.run_snakefile(args, config=config)

    if not status:
        cerr('[WARNING: full run of discovery variant calling '
             'did not successfully complete]')
    cerr(f'[Finish full run of discovery variant calling (time: {elapsed_time})]')


def main(args):
    run_multistep_variant_valler(args)

# EOF
