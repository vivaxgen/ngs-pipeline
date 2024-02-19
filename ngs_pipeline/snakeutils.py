__copyright__ = '''
snakeutils.py - ngs-pipeline
[https://github.com/vivaxgen/ngs-pipeline]

(c) 2023 Hidayat Trimarsanto <trimarsanto@gmail.com>

All right reserved.
This software is licensed under MIT license.
Please read the README.txt of this software.
'''

import os
import pathlib
import datetime
from ngs_pipeline import (cerr, cexit, arg_parser,
                          check_NGSENV_BASEDIR, check_NGS_PIPELINE_BASE,
                          get_snakefile_path, setup_config)


def init_argparser(desc, p=None):
    """ provide common arguments for snakemake-based cli"""
    p = p if p else arg_parser(desc=desc)
    p.arg_dict = {}

    # snakemake arguments
    p.add_argument('-j', type=int, default=32,
                   help='number of jobs to be executed in parallel')
    p.add_argument('--dryrun', default=False, action='store_true')
    p.add_argument('--showcmds', default=False, action='store_true')
    p.add_argument('--showconfigfiles', default=False, action='store_true')
    p.add_argument('--unlock', default=False, action='store_true')
    p.add_argument('--rerun', default=False, action='store_true')
    p.add_argument('--touch', default=False, action='store_true',
                   help='touch all output files, to avoid re-running the ngs-pipeline '
                   'such as after modifying/debugging snakemake file')
    p.add_argument('--profile', default=None,
                   help='snakemake profile to be used, also can be set from SNAKEMAKE_PROFILE env')
    p.add_argument('--nocluster', default=False, action='store_true',
                   help='run without cluster support (eg. only on local node), useful for debugging')

    # general options
    p.arg_dict['target'] = p.add_argument('--target', default='all',
                                          help='target rule in the snakefile [all]')
    p.arg_dict['snakefile'] = p.add_argument('--snakefile', default=None,
                                             help='snakemake file to be called')
    p.add_argument('-c', '--config', default=[], action='append',
                   help='config file(s) to append')
    p.add_argument('-f', '--force', default=False, action='store_true',
                   help='force the processing even if the working directory is not '
                        'under current pipeline environment base directory')
    p.add_argument('--no-config-cascade', default=False, action='store_true',
                   help='prevent from reading cascading configuration file')

    return p


def run_snakefile(args, config = {}):
    """ exeute snakefile based on args"""

    import snakemake
    import yaml

    # check sanity
    if not args.snakefile:
        cexit('ERR: Please provide snakefile to ecexute using --snakefile argument.')

    NGS_PIPELINE_BASE = check_NGS_PIPELINE_BASE()
    NGSENV_BASEDIR = pathlib.Path(check_NGSENV_BASEDIR())

    # check sanity

    cwd = pathlib.Path.cwd()
    if not args.force and not cwd.is_relative_to(NGSENV_BASEDIR):
        cexit(f'ERROR: current directory {cwd} is not relative to {NGSENV_BASEDIR}')

    configfiles = list(reversed(args.config))

    if args.no_config_cascade:
        configfiles.append(NGSENV_BASEDIR / 'config.yaml')
    else:
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

    if args.showconfigfiles:
        cerr('Config files to read:')
        cerr(yaml.dump(configfiles))
        cexit('\n')


    # for profile purposes

    if 'SNAKEMAKE_PROFILE' in os.environ:
        if args.profile is None:
            args.profile = os.environ['SNAKEMAKE_PROFILE']

    if args.profile is not None and not args.nocluster:

        profiles = [args.profile]
        parser = snakemake.get_argument_parser(profiles=profiles)
        profile_args = parser.parse_args(['--profile', args.profile])

        for k in ['cluster', 'cluster_status', 'cluster_cancel', 'jobscript']:
            arg = getattr(profile_args, k)
            if arg:
                path = pathlib.Path(arg)
                if path.exists() or path.is_absolute():
                    setattr(args, k, arg)
                else:
                    # adjust file path
                    setattr(args, k, snakemake.get_profile_file(args.profile, arg, return_default=True))
        # print(profile_args)
        # print(args)
    
    else:
        for k in ['cluster', 'cluster_status', 'cluster_cancel', 'jobscript']:
            setattr(args, k, None)
    
    # run snakefile
    start_time = datetime.datetime.now()
    status = snakemake.snakemake(
        get_snakefile_path(args.snakefile, pathlib.Path(NGS_PIPELINE_BASE) / 'rules'),
        configfiles=configfiles,
        config=setup_config(config),
        dryrun=args.dryrun,
        printshellcmds=args.showcmds,
        unlock=args.unlock,
        force_incomplete=args.rerun,
        touch=args.touch,
        cores=args.j,
        cluster=(args.cluster or os.environ.get('JOBCMD', '')) if not args.nocluster else None,
        cluster_status=args.cluster_status,
        cluster_cancel=args.cluster_cancel,
        jobscript=args.jobscript,
        targets=[args.target],
    )
    finish_time = datetime.datetime.now()

    return (status, finish_time - start_time)

# EOF

    
