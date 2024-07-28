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
    p.add_argument('-P', default=None,
                   help='procfile containing number of jobs/task')
    p.add_argument('--joblog', default=None,
                   help='name of job log file, if this is an existing directory, '
                   'the directory name will be prepended to default log file name '
                   ' [$USER-run-DATETIME.log]')
    p.add_argument('--dryrun', default=False, action='store_true')
    p.add_argument('--count', type=int, default=-1,
                   help='number of samples to be processed, useful to check initial run [-1]')
    p.add_argument('--showcmds', default=False, action='store_true')
    p.add_argument('--showconfigfiles', default=False, action='store_true')
    p.add_argument('--unlock', default=False, action='store_true')
    p.add_argument('--rerun', default=False, action='store_true')
    p.add_argument('--touch', default=False, action='store_true',
                   help='touch all output files, to avoid re-running the ngs-pipeline '
                   'such as after modifying/debugging snakemake file')
    p.add_argument('--force', default=False, action='store_true',
                   help='force running the process outside the NGSENV_BASEDIR')
    p.add_argument('--no-config-cascade', default=False, action='store_true',
                   help='prevent from doing cascading configuration')
    p.add_argument('--snakefile', default=None, choices=['var_call.smk', 'var_call_ont.smk'],
                   help='snakemake file to be run (or from VARCALL_SMK env) [var_call.smk]')
    p.add_argument('--target', choices=['all', 'mapping', 'clean'],
                   default='all',
                   help="target of snakemake module, use 'all' for GATK joint "
                   "variant or 'mapping' for FreeBayes joint variant call")
    p.add_argument('indirs', nargs='+',
                   help="directory(ies) containing the sample directory")
    return p


def run_sample_variant_caller(args):

    import pathlib
    import subprocess
    import time
    import datetime
    import getpass

    # check for number of jobs from environment
    if args.j < 0:
        if (JOBS := os.environ.get('JOBS', -1)):
            args.j = int(JOBS)
        else:
            args.j = 32

    if args.P:
        args.j = args.P

    # get NGS_PIPELINE_BASE
    NGS_PIPELINE_BASE = os.environ["NGS_PIPELINE_BASE"]

    # get NGSENV_BASEDIR
    NGSENV_BASEDIR = os.environ["NGSENV_BASEDIR"]

    # get snakefile to run
    if args.snakefile is None:
        if 'VARCALL_SMK' in os.environ:
            args.snakefile = os.environ['VARCALL_SMK']
        else:
            args.snakefile = 'var_call.smk'
    cerr(f'Snakefile to be run: {args.snakefile}')

    # set up joblog filename
    default_joblog = getpass.getuser() + '-run-' + time.strftime("%y%m%d-%H%M") + '.log'
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
        if not (args.force or indir.resolve().is_relative_to(NGSENV_BASEDIR)):
            cexit(f'[ERROR: {indir} is not relative to {NGSENV_BASEDIR}]')

        curr_samples = [s.as_posix() for s in indir.iterdir() if s.is_dir()]
        cerr(f'[Collecting {len(curr_samples)} directories from {indir}]')
        samples += curr_samples

    cerr(f'[Collecting total of {len(samples)} sample directories from '
         f'{len(args.indirs)} input directories]')

    if args.showconfigfiles:
        # only use 1 sample for showing cascading config files
        args.count = 1

    if args.count > 0:
        samples = samples[:args.count]
        cerr(f'[Limiting processing to {args.count} sample(s)]')

    # import IPython; IPython.embed()

    # sort by descending file size
    file_sizes = []
    for sample in samples:
        try:
            file_size = int(open(f'{sample}/reads/filesize').read())
        except FileNotFoundError:
            file_size = 0
        file_sizes.append(file_size)
    samples = list(zip(*(sorted(zip(file_sizes, samples), reverse=True))))[1]

    start_time = datetime.datetime.now()
    cmds = ['parallel',
            '--termseq', 'INT,1000,KILL,25',
            '--eta',
            '-j', str(args.j),
            '--joblog', args.joblog,
            '--workdir',
            '{}',
            f'ngs-pl run-indv-varcall -j 48 '
            f'{"--unlock" if args.unlock else ""} '
            f'{"--rerun" if args.rerun else ""} '
            f'{"--touch" if args.touch else ""} '
            f'{"--showcmds" if args.showcmds else ""} '
            f'{"--showconfigfiles" if args.showconfigfiles else ""} '
            f'{"--force" if args.force else ""} '
            f'{"--no-config-cascade" if args.no_config_cascade else ""} '
            f'--snakefile {args.snakefile} '
            f'{args.target}',
            ':::'
            ]

    # append sample directory list
    cmds += samples

    # run command and wait
    subprocess.call(cmds)

    cerr('\n====================== SAMPLE VARIANT CALLING RUN REPORT ======================\n')
    finish_time = datetime.datetime.now()
    cerr(f'[run-sample-variant-caller was running for {finish_time - start_time}]')

    cwd = pathlib.Path.cwd()

    # iterating the directory list and check for existance of '.finished' or '.failed'
    completed = 0
    uncompleted = []
    unknown = []
    for sample_dir in samples:
        sample_path = pathlib.Path(sample_dir)
        if (sample_path / '.completed').is_file():
            completed += 1
        elif (sample_path / '.uncompleted').is_file():
            uncompleted.append(sample_dir)
        else:
            unknown.append(sample_dir)

    cerr('')
    if any(uncompleted):
        cerr('\n  - '.join(['Uncompleted:'] + uncompleted))
        cerr('')
    if any(unknown):
        cerr('\n  - '.join(['Unknown:'] + unknown))
        cerr('')

    cerr(f'[Completed: {completed}, Uncompleted: {len(uncompleted)}, Unknown: {len(unknown)}]')
    cerr('\n')


def main(args):
    run_sample_variant_caller(args)


# EOF
