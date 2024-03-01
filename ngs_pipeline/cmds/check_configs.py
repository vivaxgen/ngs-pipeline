# check_configs.py - ngs-pipeline command
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

# to improve the responsiveness during bash autocomplete, do not import heavy
# modules (such as numpy, pandas, etc) here, but instead import them within the
# functions that require the respective heavy modules


import os
from ngs_pipeline import (cerr, cexit, arg_parser,
                          check_NGSENV_BASEDIR, check_NGS_PIPELINE_BASE,
                          get_snakefile_path, setup_config, snakeutils)


def init_argparser():
    p = snakeutils.init_argparser(desc='check configuration')
    p.arg_dict['snakefile'].choices = [
        'check_configs.smk',
    ]

    # input/output options
    p.add_argument('-o', '--outfile', required=True,
                   help='output file name where YAML will be written')

    p.add_argument('indir',
                   help='Target directory to start searching for configuration')

    return p


def check_configs(args):

    import pathlib

    args.snakefile = 'check_configs.smk'

    config = dict(outfile=args.outfile)
    status, elapsed_time = snakeutils.run_snakefile(
        args,
        config=config,
        workdir=pathlib.Path(args.indir).absolute(),
        show_configfiles=True
    )

    if not status:
        cerr('[WARNING: cannot perform configuration check properly]')
    else:
        cerr(f'[Configuration written as YAML file at {args.outfile}]')


def main(args):
    check_configs(args)

# EOF
