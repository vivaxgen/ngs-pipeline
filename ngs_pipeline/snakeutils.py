__copyright__ = '''
snakeutils.py - ngs-pipeline
[https://github.com/vivaxgen/ngs-pipeline]

(c) 2023-2024 Hidayat Trimarsanto <trimarsanto@gmail.com>

All right reserved.
This software is licensed under MIT license.
Please read the README.txt of this software.
'''

import os
import sys
import pathlib
import datetime
import ngs_pipeline
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

    # configuration files
    p.add_argument('--base-config', default=None,
                   help='path for base configuration file, relative to '
                   'base environment directory')
    p.add_argument('-p', '--panel', default=None,
                   help='panel to be used (eg. PANEL -> configs/PANEL.yaml as base config)')
    p.add_argument('-c', '--config', default=[], action='append',
                   help='config file(s) to append')
    p.add_argument('-f', '--force', default=False, action='store_true',
                   help='force the processing even if the working directory is not '
                        'under current pipeline environment base directory')
    p.add_argument('--no-config-cascade', default=False, action='store_true',
                   help='prevent from reading cascading configuration file')

    return p


def run_snakefile(args, config={}, workdir=None,
                  show_configfiles=False):
    """ execute snakefile based on args """

    import snakemake
    import yaml

    if snakemake.__version__ >= '8':
        return run_snakefile_8(args, config, workdir, show_configfiles)

    # check sanity
    if not args.snakefile:
        cexit('ERR: Please provide snakefile to ecexute using --snakefile argument.')

    # getting values from environment
    NGS_PIPELINE_BASE = check_NGS_PIPELINE_BASE()
    NGSENV_BASEDIR = pathlib.Path(check_NGSENV_BASEDIR())
    if 'NGS_PIPELINE_FORCE' in os.environ:
        args.force = True
    if 'NGS_PIPELINE_NO_CONFIG_CASCADE' in os.environ:
        args.no_config_cascade = True

    # check sanity

    cwd = workdir or pathlib.Path.cwd()
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

    # get panel configuration and set as base configuration from configs/
    if args.panel:
        args.base_config = 'configs/' + args.panel + '.yaml'

    if args.base_config:
        configfiles.append(NGSENV_BASEDIR / args.base_config)

    if not any(configfiles):
        cexit(f'ERROR: cannot find any config.yaml in {config_dirs}')

    configfiles.reverse()

    if args.showconfigfiles:
        cerr('Config files to read:')
        cerr(yaml.dump([cf.as_posix() for cf in configfiles]))
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
                    setattr(args, k,
                            snakemake.get_profile_file(args.profile,
                                                       arg,
                                                       return_default=True))
        # print(profile_args)
        # print(args)
    
    else:
        for k in ['cluster', 'cluster_status', 'cluster_cancel', 'jobscript']:
            setattr(args, k, None)
    
    # run snakefile
    start_time = datetime.datetime.now()
    status = snakemake.snakemake(
        get_snakefile_path(
            args.snakefile,
            pathlib.Path(NGS_PIPELINE_BASE) / 'ngs_pipeline' / 'rules'),
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

    if show_configfiles:
        cerr('Config files read are:')
        cerr(yaml.dump([cf.as_posix() for cf in configfiles]))

    return (status, finish_time - start_time)


def run_snakefile_8(args, config={}, workdir=None,
                  show_configfiles=False):
    """ execute snakefile based on args """

    import snakemake
    from snakemake import cli
    import yaml

    # monkey-patch snakemake
    cli_parse_config = cli.parse_config

    def parse_config(entries):
        if type(entries) == dict:
            return entries
        return cli_parse_config(entries)

    cli.parse_config = parse_config
    # end of monkey patching

    # check sanity
    if not args.snakefile:
        cexit('ERR: Please provide snakefile to ecexute using --snakefile argument.')

    # getting values from environment
    NGS_PIPELINE_BASE = check_NGS_PIPELINE_BASE()
    NGSENV_BASEDIR = pathlib.Path(check_NGSENV_BASEDIR())
    if 'NGS_PIPELINE_FORCE' in os.environ:
        args.force = True
    if 'NGS_PIPELINE_NO_CONFIG_CASCADE' in os.environ:
        args.no_config_cascade = True

    # check sanity

    cwd = workdir or pathlib.Path.cwd()
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

    # get panel configuration and set as base configuration from configs/
    if args.panel:
        args.base_config = 'configs/' + args.panel + '.yaml'

    if args.base_config:
        configfiles.append(NGSENV_BASEDIR / args.base_config)

    if not any(configfiles):
        cexit(f'ERROR: cannot find any config.yaml in {config_dirs}')

    configfiles.reverse()

    if args.showconfigfiles:
        cerr('Config files to read:')
        cerr(yaml.dump([cf.as_posix() for cf in configfiles]))
        cexit('\n')

    # for profile purposes

    if 'SNAKEMAKE_PROFILE' in os.environ:
        if args.profile is None:
            args.profile = os.environ['SNAKEMAKE_PROFILE']

    if args.nocluster:
        # nocluster means prevent from running using batch/job scheduler
        args.profile = None

    try:
        # set with --profile first
        if args.profile:
            argv = ['--profile', args.profile]
        else:
            argv = []

        print(argv)

        # XXX: need to modify to use snakemake API
        parser, cargs = cli.parse_args(argv)

        # set cargs further from args
        cargs.snakefile = get_snakefile_path(
            args.snakefile,
            from_module=ngs_pipeline
        )
        cargs.configfile = configfiles
        #cargs.config = [f'{k}={v}' for k, v in setup_config(config).items()]
        cargs.config = setup_config(config)
        cargs.targets = [args.target]

        # running mode
        cargs.dryrun = args.dryrun
        cargs.touch = args.touch
        cargs.rerun_incomplete = args.rerun
        cargs.unlock = args.unlock

        # running parameters
        cargs.cores = args.j
        cargs.printshellcmds = args.showcmds

        start_time = datetime.datetime.now()
        status = cli.args_to_api(cargs, parser)
        finish_time = datetime.datetime.now()

        return (status, finish_time - start_time)

    except Exception as e:
        cli.print_exception(e)
        sys.exit(1)

    raise RuntimeError('FATAL ERROR: should not execute this part of code')


def scan_for_config_keywords(path):
    """ return a list of keywords used as config keys in any of the
        snakemake rules
    """

    import re

    mo_bracket = re.compile(r'''config\s*\[\s*['"]([^]]*)['"]\s*\]''')
    mo_get = re.compile(r'''config\s*\.\s*get\s*\(\s*['"]([^'"]*)''')

    keywords = []

    path = pathlib.Path(path)

    for rule_file in path.glob('*.smk'):
        with open(rule_file) as f_in:
            for line in f_in:
                keys = mo_bracket.findall(line) + mo_get.findall(line)
                if any(keys):
                    keywords += keys

    return sorted(set(keywords))


# EOF
