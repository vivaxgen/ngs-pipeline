#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

__copyright__ = '''
run-targeted-varcall.py - ngs-pipeline command line
[https://github.com/vivaxgen/ngs-pipeline]

(c) 2023 Hidayat Trimarsanto <trimarsanto@gmail.com>

All right reserved.
This software is licensed under MIT license.
Please read the README.txt of this software.
'''

# to improve the responsiveness during bash autocomplete, do not import heavy
# modules (such as numpy, pandas, etc) here, but instead import them within the
# functions that requires respective heavy modules

import sys
import os

# check that we have NGS_PIPELINE_BASE environemt
if 'NGS_PIPELINE_BASE' not in os.environ:
    print('ERROR: NGS_PIPELINE_BASE environment is not set. '
          'Please set proper shell enviroment by sourcing relevant activate.sh',
          file=sys.stderr)
    sys.exit(1)

# check that we have NGSENV_BASEDIR environemt
if 'NGSENV_BASEDIR' not in os.environ:
    print('ERROR: NGSENV_BASEDIR environment is not set. '
          'Please set proper shell enviroment by sourcing relevant activate.sh',
          file=sys.stderr)
    sys.exit(1)


from ngsutils import cerr, cexit, run_main, arg_parser


def init_argparser():
    p = arg_parser(desc='run targeted variant calling')
    p.add_argument('-j', type=int, default=32)

    # Snakemake options
    p.add_argument('--dryrun', default=False, action='store_true')
    p.add_argument('--showcmds', default=False, action='store_true')
    p.add_argument('--unlock', default=False, action='store_true')
    p.add_argument('--rerun', default=False, action='store_true')

    # general options
    p.add_argument('--target', default='all',
                   choices=['all'],
                   help='target rule in the snakefile [all]')
    p.add_argument('-o', '--outdir', default='analysis',
                   help='directory for output [analysis/]')
    p.add_argument('--snakefile', default='targeted_varcall_pe.smk',
                   choices=['targeted_varcall_pe.smk', 'targeted_varcall_lr.smk'],
                   help='snakemake file to be called [targeted_varcall_pe.smk]')
    p.add_argument('-c', '--config', default=[], action='append',
                   help='config file(s) to append')

    # input options
    p.add_argument('infiles', nargs='+',
                   help='FASTQ input files, eg. sample-1.fastq.gz')

    return p


def run_targeted_varcall(args):

    import pathlib
    import snakemake
    import datetime

    cwd = pathlib.Path.cwd()
    NGSENV_BASEDIR = pathlib.Path(os.environ['NGSENV_BASEDIR'])

    # check sanity
    if not cwd.is_relative_to(NGSENV_BASEDIR):
        cexit(f'ERROR: current directory {cwd} is not relative to {NGSENV_BASEDIR}')

    configfiles = []

    # for each config directory, check config file existence
    config_dirs = []
    config_path = cwd
    while config_path.is_relative_to(NGSENV_BASEDIR):
        config_dirs.append(config_path)
        configfile = config_path / 'config.yaml'
        if configfile.is_file():
            configfiles.append(configfile)
        config_path = config_path.parent

    if not any(configfiles):
        cexit(f'ERROR: cannot find any config.yaml in {config_dirs}')

    configfiles.reverse()

    # run snakemake

    start_time = datetime.datetime.now()
    status = snakemake.snakemake(
        pathlib.Path(os.environ['NGS_PIPELINE_BASE']) / 'smk' / args.snakefile,
        configfiles=configfiles,
        config=dict(infiles=args.infiles,
                    outdir=args.outdir),
        dryrun=args.dryrun,
        printshellcmds=args.showcmds,
        unlock=args.unlock,
        force_incomplete=args.rerun,
        cores=args.j,
        cluster=os.environ.get('JOBCMD', ''),
        cluster_cancel="scancel",
        targets=[args.target],
    )
    finish_time = datetime.datetime.now()
    if not status:
        cerr('[WARNING: targeted variant calling did not successfully complete]')
    cerr(f'[Finish targeted variant calling (time: {finish_time - start_time})]')


if __name__ == '__main__':
    run_main(init_argparser, run_targeted_varcall)

# EOF
