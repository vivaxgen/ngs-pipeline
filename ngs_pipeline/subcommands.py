# subcommands.py
# [https://github.com/trmznt/subcommands]

__copyright__ = "(c) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"
__version__ = '2024.03.04'

# this module provides subcommands, eg. PROG subcommand [options]

import os
import sys
import platform
import pathlib
import importlib
import argparse
import argcomplete
from typing import Callable


def _cout(msg: str):
    print(msg, file=sys.stdout)


def _cerr(msg: str):
    print(msg, file=sys.stderr)
    sys.stderr.flush()


def _cexit(msg: str, exit_code: int = 1):
    _cerr(msg)
    sys.exit(exit_code)


def arg_parser(desc: str = ''):
    p = argparse.ArgumentParser(description=desc)
    return add_debug_to_parser(p)


def add_debug_to_parser(p: argparse.ArgumentParser):

    try:
        p.add_argument('--debug', action='store_true', default=False,
                       help='open ipdb when uncatched exception is raised')
    except argparse.ArgumentError:
        pass

    return p


class SubCommands(object):

    def __init__(
        self,
        module_env: str | None = None,
        modules: list[str] = [],
        allow_any_script: bool = False,
        allow_shell: bool = False,
        greet_func: Callable | None = None,
        usage_func: Callable | None = None,
        help_func: Callable | None = None
    ):

        self.allow_any_script = allow_any_script
        self.allow_shell = allow_shell
        self.greet = greet_func or self.generic_greet
        self.usage = usage_func or self.generic_usage
        self.help = help_func or self.generic_help
        self.prog_name = sys.argv[0].rsplit('/', 1)[-1]

        # set up module list
        self.modules = modules
        if module_env:
            self.modules += [m for m in os.environ[module_env].split(':') if m]

    def autocomplete(self, tokens: list[str]):

        last_token = tokens[-1]

        # prepare the completion lists
        cmds = self.get_command_list()

        if len(tokens) == 1:
            completions = sorted(cmds)

        else:
            completions = sorted(
                [opt for opt in cmds if opt.startswith(last_token)]
            )

        # send results through fd 8
        ifs = os.environ.get('IFS', '\013')
        out_stream = os.fdopen(8, 'w')
        out_stream.write(ifs.join(completions))
        sys.exit(0)

    def generic_usage(self):
        _cexit(
            f'  usage:\n'
            f'      {self.prog_name} CMD [OPTIONS]\n'
            f'      {self.prog_name} [-i] [-l] [-h]\n'
            f'  try: {self.prog_name} -l'
        )

    def generic_greet(self):
        _cerr(f'{self.prog_name} - command line interface')
        _cerr(f'Host: {platform.uname().node}')

    def generic_help(self):
        self.generic_usage()

    def get_command_list(self) -> list[str]:
        # read sqpy.cmds directory
        cmd_files = []
        for module in self.modules:
            M = importlib.import_module(module)
            cmd_directory = pathlib.Path(M.__path__[0])
            cmd_files += cmd_directory.iterdir()
        cmds = set(
            [p.name.removesuffix('.py').replace('_', '-') for p in cmd_files]
        ) - {'--init--', '--pycache--'}
        return sorted(cmds)

    def show_commands(self):
        _cout('Available commands:')
        for cmd in self.get_command_list():
            _cout(f'  {cmd}')
        sys.exit(0)

    def run_cmd(self, args: list[str]):
        """ this method will run module init_argparser() and main() """

        command = args[0].replace('-', '_')
        for module in self.modules:
            try:
                module_name = f'{module}.{command}'
                M = importlib.import_module(module_name)
                break
            except ModuleNotFoundError as err:
                if err.name != module_name:
                    raise
        else:
            _cexit(f'Command: {args[0]} does not exist! '
                   f'Please check with "{self.prog_name} -l".')

        _cerr(f'Executing: {M.__file__}')
        parser = M.init_argparser()
        if not parser:
            _cexit('Fatal ERR: init_argparser() does not return properly')
        add_debug_to_parser(parser)

        # perform bash autocompletion if needed
        argcomplete.autocomplete(parser)

        # if autocomplete does not exit:
        args = parser.parse_args(args[1:])
        if args.debug:
            from ipdb import launch_ipdb_on_exception
            with launch_ipdb_on_exception():
                _cerr('WARN: running in debug mode')
                M.main(args)
        else:
            M.main(args)

    def run_script(self, args: list[str]):
        """ this method will run an arbitrary python script """

        path = args[0]
        _cerr(f'Attempting to run script: {path}')

        # expand home directory
        if path.startswith('~'):
            path = pathlib.Path(path).expanduser()

        with open(path) as fh:
            code = compile(fh.read(), path, 'exec')
            _l = {'__name__': '__anyscript_main__'}
            exec(code, None, _l)
            if 'main' in _l:
                globals().update(_l)
                main = _l['main']
                if 'init_argparser' in _l:
                    init_argparser = _l['init_argparser']
                    p = init_argparser()
                    if not isinstance(p, argparse.ArgumentParser):
                        _cexit(
                            'ERR: init_argparser() did not return '
                            'ArgumentParser instance')
                    add_debug_to_parser(p)

                    argcomplete.autocomplete(p)
                    argp = p.parse_args(args)

                    if argp.debug:
                        from ipdb import launch_ipdb_on_exception
                        with launch_ipdb_on_exception():
                            _cerr('WARN: running in debug mode')
                            main(argp)
                    else:
                        main(argp)
                else:
                    import inspect
                    if 'args' in inspect.signature(main).parameters:
                        main(args=args)
                    else:
                        main()

    def main(self):

        tokens = []

        if '_ARGCOMPLETE' in os.environ:
            # we are invoked from argcomplete
            line = os.environ.get('COMP_LINE', '')
            tokens = line.split()
            if len(tokens) == 1 or (len(tokens) == 2 and
                                    not line.endswith(' ')):
                self.autocomplete(tokens)
            os.environ['COMP_LINE'] = line.split(' ', 1)[1]
            os.environ['COMP_POINT'] = str(len(os.environ['COMP_LINE']))
            cmd = tokens[1]
        else:
            self.greet()
            if len(sys.argv) == 1:
                self.usage()
            cmd = sys.argv[1]

        if cmd == '-i' and self.allow_shell:
            # go to interactive mode
            import IPython
            IPython.embed()

        if cmd == '-l':
            # show available commands
            self.show_commands()

        elif cmd == '-h':
            # show some help
            self.help()

        elif cmd.endswith('.py') and self.allow_any_script:

            # run script, with:
            #   path = cmd
            #   args = sys.argv[2:] (eg. 'my_prog my_script arg0 arg1')
            self.run_script(cmd, sys.argv[2:])

        else:
            self.run_cmd(tokens[1:] if any(tokens) else sys.argv[1:])


# EOF
