#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

__copyright__ = '''
run_varcall.py - ngs-pipeline command line
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

from ngs_pipeline import cerr, cexit, run_main, arg_parser


# usage: run_varcall.py

def init_argparser():
    p = arg_parser(desc='run varcalling pipeline')
    p.add_argument('-j', type=int, default=72)
    p.add_argument('--dryrun', default=False, action='store_true')
    p.add_argument('--showcmds', default=False, action='store_true')
    p.add_argument('--unlock', default=False, action='store_true')
    p.add_argument('--rerun', default=False, action='store_true')
    p.add_argument('--touch', default=False, action='store_true',
                   help='touch all output files, to avoid re-running the ngs-pipeline '
                   'such as after modifying/debugging snakemake file')
    p.add_argument('--snakefile', default='var_call.smk',
                   help='snakemake file to be run [var_call.smk]')
    p.add_argument('target')
    p.add_argument('sampledir')
    return p


def run_varcall(args):

    import pathlib
    import snakemake
    import datetime

    # directory structure:
    #  [NGSENV_BASEDIR]/sets/[BATCH]/samples/[SAMPLE]
    #
    # possible configfiles:
    #  [NGSENV_BASEDIR]/sets/[BATCH]/samples/[SAMPLE]/config.yaml
    #  [NGSENV_BASEDIR]/sets/[BATCH]/samples/config.yaml
    #  [NGSENV_BASEDIR]/sets/[BATCH]/config.yaml
    #  [NGSENV_BASEDIR]/config.yaml

    cwd = pathlib.Path.cwd()
    NGSENV_BASEDIR = pathlib.Path(os.environ['NGSENV_BASEDIR'])

    # check sanity
    if not cwd.is_relative_to(NGSENV_BASEDIR):
        cexit(f'ERROR: current directory {cwd} is not relative to {NGSENV_BASEDIR}')

    # sanity check
    #if (cwd_part := cwd.parts[-1]) != args.sample:
    #    cexit(f'ERROR: sample argument {args.sample} is not identical with '
    #          f'last part of current directory name {cwd_part}!')
    sample = cwd.parts[-1]

    if not (cwd / 'reads').is_dir():
        cexit(f'ERROR: current directory {cwd} does not have "reads" directory!')

    configfiles = []

    # for each config directory, check config file existence
    #config_dirs = [cwd, cwd.parent, cwd.parent.parent, cwd.parent.parent.parent.parent]
    config_dirs = []
    config_path = cwd
    #for config_path in config_dirs:
    while config_path.is_relative_to(NGSENV_BASEDIR):
        config_dirs.append(config_path)
        configfile = config_path / 'config.yaml'
        if configfile.is_file():
            configfiles.append(configfile)
        config_path = config_path.parent

    if not any(configfiles):
        cexit(f'ERROR: cannot find any config.yaml in {config_dirs}')

    configfiles.reverse()

    # remove .failed or .finished screen first if exists
    (cwd / '.completed').unlink(missing_ok=True)
    (cwd / '.uncompleted').unlink(missing_ok=True)

    # run smk

    start_time = datetime.datetime.now()
    status = snakemake.snakemake(
        pathlib.Path(os.environ['NGS_PIPELINE_BASE']) / 'rules' / args.snakefile,
        configfiles=configfiles,
        dryrun=args.dryrun,
        touch=args.touch,
        printshellcmds=args.showcmds,
        unlock=args.unlock,
        force_incomplete=args.rerun,
        cores=args.j,
        cluster=os.environ.get('JOBCMD', ''),
#        cluster="sbatch --cpus-per-task={threads} --job-name=smk-{rule}-{sample}",
        cluster_cancel="scancel",
        targets=[args.target],
    )
    finish_time = datetime.datetime.now()
    open(cwd / ('.uncompleted' if not status else '.completed'), 'w').write(f'{finish_time - start_time}')


if __name__ == '__main__':
    run_main(init_argparser, run_varcall)

# EOF
