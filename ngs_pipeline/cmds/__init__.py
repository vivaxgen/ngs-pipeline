__copyright__ = '''
ngs_pipeline/cmds/__init__.py
[https://github.com/vivaxgen/ngs-pipeline]

(c) 2023 Hidayat Trimarsanto <trimarsanto@gmail.com>

All right reserved.
This software is licensed under MIT license.
Please read the README.txt of this software.
'''


from ngs_pipeline import cerr, cexit, get_command_modules
import argcomplete
import importlib
import pathlib
import os


def run_main(args):

    command = args[0].replace('-', '_')
    for module in get_command_modules():
        try:
            module_name = f'{module}.{command}'
            M = importlib.import_module(module_name)
            break
        except ModuleNotFoundError as err:
            if err.name != module_name:
                raise
    else:
        cexit(f'Command: {args[0]} does not exist! Please check with "ngs-pl showcmds".')

    cerr(f'Executing: {M.__file__}')
    parser = M.init_argparser()
    if not parser:
        cexit('Fatal ERR: init_argparser() does not return properly')

    # perform bash autocompletion if needed
    argcomplete.autocomplete(parser)

    # if autocomplete doe not exit:
    args = parser.parse_args(args[1:])
    if args.debug:
        from ipdb import launch_ipdb_on_exception
        with launch_ipdb_on_exception():
            cerr('WARN: running in debug mode')
            M.main(args)
    else:
        M.main(args)


def list_commands():
    # read sqpy.cmds directory
    cmd_files = []
    for  module in get_command_modules():
        M = importlib.import_module(module)
        cmd_directory = pathlib.Path(M.__path__[0])
        cmd_files += cmd_directory.iterdir()
    cmds = set(
        [p.name.removesuffix('.py').replace('_', '-') for p in cmd_files]
    ) - {'--init--', '--pycache--'}
    return sorted(cmds)

# EOF
