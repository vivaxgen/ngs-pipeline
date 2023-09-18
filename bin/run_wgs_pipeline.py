#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

__copyright__ = '''
run_wgs_pipeline.py - ngs-pipeline command line
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
    print('ERROR: please set proper shell enviroment by sourcing activate.sh',
          file=sys.stderr)
    sys.exit(1)


from ngsutils import cerr, cexit, run_main, arg_parser


# usage: run_varcall.py

def init_argparser():
    p = arg_parser(desc='run WGS mapping & genotyping pipeline')
    p.add_argument('-j', type=int, default=-1,
                   help='number of samples to be processed in parallel, will override '
                   'JOBS environment [32]')
    p.add_argument('--joblog', default=None,
                   help='name of job log file [run-DATETIME.log]')
    p.add_argument('--dryrun', default=False, action='store_true')
    p.add_argument('--count', type=int, default=-1,
                   help='number of samples to be processed, useful to check initial run [-1]')
    p.add_argument('--showcmds', default=False, action='store_true')
    p.add_argument('--unlock', default=False, action='store_true')
    p.add_argument('--rerun', default=False, action='store_true')
    p.add_argument('--touch', default=False, action='store_true',
                   help='touch all output files, to avoid re-running the ngs-pipeline '
                   'such as after modifying/debugging snakemake file')
    p.add_argument('--target', choices=['all', 'mapping', 'clean'],
                   help="target of snakemake module, use 'all' for GATK joint variant "
                   "or 'mapping' for FreeBayes joint variant call")
    p.add_argument('indirs', nargs='+',
                   help="directory(ies) containing the sample directory")
    return p


def run_wgs_pipeline(args):

    import pathlib
    import subprocess

    # check for number of jobs from environment
    if args.j < 0:
        if (JOBS := os.environ.get('JOBS', -1)):
            args.j = int(JOBS)
        else:
            args.j = 32

    # get NGS_PIPELINE_BASE
    NGS_PIPELINE_BASE = os.environ["NGS_PIPELINE_BASE"]

    # get NGSENV_BASEDIR
    NGSENV_BASEDIR = os.environ["NGSENV_BASEDIR"]

    # joblog file
    if not args.joblog:
        args.joblog = 'run-' + time.strftime("%y%m%d-%H%M") + '.log'

    cwd = pathlib.Path.cwd()

    # look for sample directories

    samples = []
    for indir in args.indirs:
        indir = pathlib.Path(indir).resolve()

        #import IPython; IPython.embed()

        # check whether indir is a directory
        if not indir.is_dir():
            cexit(f'[ERROR: directory {indir} does not exist]')

        # check whether indir is under NGS_PIPELINE_BASE
        if not indir.resolve().is_relative_to(NGSENV_BASEDIR):
            cexit(f'[ERROR: {indir} is not relative to {NGSENV_BASEDIR}]')

        curr_samples = [s.as_posix() for s in indir.iterdir() if s.is_dir()]
        cerr(f'[Collecting {len(curr_samples)} directories from {indir}]')
        samples += curr_samples
    cerr(f'[Collecting total of {len(samples)} sample directories from '
         f'{len(args.indirs)} input directories]')

    if args.count > 0:
        samples = samples[:args.count]
        cerr(f'[Limiting processing to {args.count} sample(s)]')

    # import IPython; IPython.embed()

    cmds = ['parallel',
            '--eta',
            '-j', str(args.j),
            '--joblog', args.joblog,
            '--workdir',
            '{}',
            f'{NGS_PIPELINE_BASE}/bin/run_varcall.py -j 48 '
            f'{"--unlock" if args.unlock else ""} '
            f'{"--rerun" if args.rerun else ""} '
            f'{"--touch" if args.touch else ""} '
            f'{args.target}',
            ':::'
            ]

    # append sample directory list

    cmds += samples

    subprocess.call(cmds)

    cwd = pathlib.Path.cwd()


if __name__ == '__main__':
    run_main(init_argparser, run_wgs_pipeline)

# EOF
