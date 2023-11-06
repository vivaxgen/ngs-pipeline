__copyright__ = '''
run_sample_variant_caller.py - ngs-pipeline command line
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
from ngs_pipeline import cerr, cexit, arg_parser


# usage: run_varcall.py

def init_argparser():
    p = arg_parser(desc='run WGS mapping & genotyping pipeline')
    p.add_argument('-j', type=int, default=-1,
                   help='number of samples to be processed in parallel, will override '
                   'JOBS environment [32]')
    p.add_argument('--joblog', default=None,
                   help='name of job log file, if this is an existing directory, '
                   'the directory name will be prepended to default log file name '
                   ' [$USER-run-DATETIME.log]')
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


def run_sample_variant_caller(args):

    import pathlib
    import subprocess
    import time
    import datetime

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

    # set up joblog filename
    default_joblog = os.getlogin() + '-run' + time.strftime("%y%m%d-%H%M") + '.log'
    if args.joblog:
        if (jobpath := pathlib.Path(args.joblog)).is_dir():
            args.joblog = jobpath / default_joblog
        else:
            args.joblog = jobpath
    else:
        args.joblog = pathlib.Path(default_joblog)

    cwd = pathlib.Path.cwd()

    # look for sample directories

    samples = []
    for indir in args.indirs:
        indir = pathlib.Path(indir).resolve()

        #import IPython; IPython.embed()

        # check whether indir is a directory
        if not indir.is_dir():
            cexit(f'[ERROR: directory {indir} does not exist]')

        # check whether indir is under NGSENV_BASEDIR
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

    start_time = datetime.datetime.now()
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

    finish_time = datetime.datetime.now()
    cerr(f'[WGS pipeline was running for {finish_time - start_time}]')

    cwd = pathlib.Path.cwd()

    # iterating the directory list and check for existance of '.finished' or '.failed'
    finished = 0
    failed = []
    unknown = []
    for sample_dir in samples:
        sample_path = pathlib.Path(sample_dir)
        if (sample_path / '.finished').is_file():
            finished += 1
        elif (sample_path / '.failed').is_file():
            failed.append(sample_dir)
        else:
            unknown.append(sample_dir)

    cerr(f'[Finished: {finished}, Failed: {len(failed)}, Unknown: {len(unknown)}]')
    if any(failed):
        cerr('\n  - '.join(['Failed:'] + failed))
    if any(unknown):
        cerr('\n  - '.join(['Unknown:'] + unknown))


def main(args):
    run_sample_variant_caller(args)


# EOF
