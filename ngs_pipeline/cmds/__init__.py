__copyright__ = '''
ngs_pipeline/cmds/__init__.py
[https://github.com/vivaxgen/ngs-pipeline]

(c) 2023 Hidayat Trimarsanto <trimarsanto@gmail.com>

All right reserved.
This software is licensed under MIT license.
Please read the README.txt of this software.
'''


from ngs_pipeline import cerr, cexit
import argcomplete
import importlib
import pathlib


def run_main(args):

    command = args[0].replace('-', '_')
    M = importlib.import_module('ngs_pipeline.cmds.' + command)
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
    import ngs_pipeline.cmds
    cmds_directory = pathlib.Path(ngs_pipeline.cmds.__file__).parent
    cmds = set(
        [p.name.removesuffix('.py').replace('_', '-') for p in cmds_directory.iterdir()]
    ) - {'--init--', '--pycache--'}
    return sorted(cmds)

# EOF
