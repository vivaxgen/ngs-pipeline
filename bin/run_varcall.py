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

from ngsutils import cerr, cexit, run_main, arg_parser


# usage: run_varcall.py

def init_argparser():
    p = arg_parser(desc='run varcalling pipeline')
    p.add_argument('-j', type=int, default=72)
    p.add_argument('--dryrun', default=False, action='store_true')
    p.add_argument('--showcmds', default=False, action='store_true')
    p.add_argument('--unlock', default=False, action='store_true')
    p.add_argument('--rerun', default=False, action='store_true')
    p.add_argument('target')
    p.add_argument('sampledir')
    return p


def run_varcall(args):

    import pathlib
    import snakemake

    # directory structure:
    #  [NGSENV_BASEDIR]/sets/[BATCH]/samples/[SAMPLE]
    #
    # possible configfiles:
    #  [NGSENV_BASEDIR]/sets/[BATCH]/samples/[SAMPLE]/config.yaml
    #  [NGSENV_BASEDIR]/sets/[BATCH]/samples/config.yaml
    #  [NGSENV_BASEDIR]/sets/[BATCH]/config.yaml
    #  [NGSENV_BASEDIR]/config.yaml

    cwd = pathlib.Path.cwd()

    # sanity check
    #if (cwd_part := cwd.parts[-1]) != args.sample:
    #    cexit(f'ERROR: sample argument {args.sample} is not identical with '
    #          f'last part of current directory name {cwd_part}!')
    sample = cwd.parts[-1]

    if not (cwd / 'reads').is_dir():
        cexit('ERROR: current directory does not have "reads" directory!')

    configfiles = []

    # for each config directory, check config file existence
    config_dirs = [cwd, cwd.parent, cwd.parent.parent, cwd.parent.parent.parent.parent]
    for config_path in config_dirs:
        configfile = config_path / 'config.yaml'
        if configfile.is_file():
            configfiles.append(configfile)

    if not any(configfiles):
        cexit(f'ERROR: cannot find any config.yaml in {config_dirs}')

    configfiles.reverse()

    # run smk

    snakemake.snakemake(
        pathlib.Path(os.environ['NGS_PIPELINE_BASE']) / 'smk' / 'var_call.smk',
        configfiles=configfiles,
        dryrun=args.dryrun,
        printshellcmds=args.showcmds,
        unlock=args.unlock,
        force_incomplete=args.rerun,
        cores=args.j,
        cluster=os.environ.get('JOBCMD', ''),
#        cluster="sbatch --cpus-per-task={threads} --job-name=smk-{rule}-{sample}",
        cluster_cancel="scancel",
        targets=[args.target],
    )


if __name__ == '__main__':
    run_main(init_argparser, run_varcall)

# EOF
