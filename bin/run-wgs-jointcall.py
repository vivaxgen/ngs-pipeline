#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

__copyright__ = '''
prepare_samples.py - ngs-pipeline command line
[https://github.com/vivaxgen/ngs-pipeline]

(c) 2022 Hidayat Trimarsanto <trimarsanto@gmail.com>

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

if 'NGSENV_BASEDIR' not in os.environ:
    print('ERROR: NGSENV_BASEDIR environment is not set. '
          'Please set proper shell enviroment by sourcing relevant activate.sh',
          file=sys.stderr)
    sys.exit(1)


from ngsutils import cerr, cexit, run_main, arg_parser


def init_argparser():
    p = arg_parser(desc='run WGS varcalling pipeline')
    p.add_argument('-j', type=int, default=72)
    p.add_argument('--dryrun', default=False, action='store_true')
    p.add_argument('--showcmds', default=False, action='store_true')
    p.add_argument('--unlock', default=False, action='store_true')
    p.add_argument('--rerun', default=False, action='store_true')
    p.add_argument('--cluster',
                   default='sbatch --cpus-per-task={threads} --job-name=smk-{rule}-{wildcards.reg} --output=.slurm-%j.out')
    p.add_argument('--nocluster', default=False, action='store_true',
                   help='run without cluster support (eg. only on local node')
    p.add_argument('--target', default='all')
    p.add_argument('source_dirs', nargs='+',
                   help='source directories containing sample directories')
    return p


def run_pipeline(args):

    import pathlib
    import snakemake
    import datetime

    # sanity check for all directory
    for srcdir in args.source_dirs:
        if not pathlib.Path(srcdir).is_dir():
            cexit(f'ERROR: directory {srcdir} does not exist!')

    cerr('[Step: joint variant calling]')
    start_time = datetime.datetime.now()
    status = snakemake.snakemake(
        pathlib.Path(os.environ['NGS_PIPELINE_BASE']) / 'smk' / 'jointvarcall_gatk.smk',
        configfiles=[pathlib.Path(os.environ['NGSENV_BASEDIR']) / 'config.yaml'],
        config=dict(srcdirs=args.source_dirs),
        printshellcmds=args.showcmds,
        dryrun=args.dryrun,
        unlock=args.unlock,
        force_incomplete=args.rerun,
        cores=args.j,
        cluster=args.cluster if not args.nocluster else None,
        cluster_cancel='scancel' if not args.nocluster else None,
        targets=[args.target],
    )
    finish_time = datetime.datetime.now()
    if not status:
        cerr('[WARNING: joint varian calling step did not successfully complete]')
    cerr(f'[Finish joint variant calling (time: {finish_time - start_time})]')


if __name__ == '__main__':
    run_main(init_argparser, run_pipeline)

# EOF

