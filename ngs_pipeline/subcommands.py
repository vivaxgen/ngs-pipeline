# subcommands.py
# [https://github.com/trmznt/subcommands]

__copyright__ = "(c) 2024-2025, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"
__version__ = "2025.10.20.01"

# this module provides subcommands, eg. PROG subcommand [options]

import os
import sys
import platform
import pathlib
import importlib
import argparse
import argcomplete
import logging
from typing import Callable

L = logging.getLogger(__name__)


def _cout(msg: str):
    print(msg, file=sys.stdout)


def _cerr(msg: str):
    print(msg, file=sys.stderr)
    sys.stderr.flush()


def _cexit(msg: str, exit_code: int = 1):
    _cerr(msg)
    sys.exit(exit_code)


__ARGUMENT_PARSER__ = argparse.ArgumentParser


def set_argument_parser_class(class_):
    global __ARGUMENT_PARSER__
    __ARGUMENT_PARSER__ = class_


def arg_parser(desc: str = ""):

    # rename program name to reflect the subcommand
    prog_name = argparse._os.path.basename(argparse._sys.argv[0])
    subprog_name = argparse._sys.argv[1]
    p = __ARGUMENT_PARSER__(prog=f"{prog_name} {subprog_name}", description=desc)
    return add_debug_to_parser(p)


def add_debug_to_parser(p: argparse.ArgumentParser):

    try:
        p.add_argument(
            "--debug",
            action="store_true",
            default=False,
            help="open ipdb when uncatched exception is raised",
        )
    except argparse.ArgumentError:
        pass

    return p


__current_subcommands__ = None


def get_subcommands():
    global __current_subcommands__
    return __current_subcommands__


class SubCommands(object):

    def __init__(
        self,
        # list of Python module names to check for commands, eg:
        # ['ngs_pipeline.cmds'] or ['seqpy.mds']
        # the order implies precedence
        modules: list[str] = [],
        # environment variable name to be used to get list of Python modules
        # to check for commands, the order implies precedence
        module_env: str | None = None,
        # whether mdules named in env variable takes precedence
        env_takes_precedence: bool = False,
        # set true to allow running scripts located in any directory
        allow_any_script: bool = False,
        # set true to allow for IPython shell
        allow_shell: bool = False,
        # set for specific greet, usage and help functions
        greet_func: Callable | None = None,
        usage_func: Callable | None = None,
        help_func: Callable | None = None,
        # additional copyright
        copyright: str | None = None,
    ):
        global __current_subcommands__
        __current_subcommands__ = self

        L.debug("initializing SubCommands class")
        self.allow_any_script = allow_any_script
        self.allow_shell = allow_shell
        self.greet = greet_func or self.generic_greet
        self.usage = usage_func or self.generic_usage
        self.help = help_func or self.generic_help
        self.prog_name = sys.argv[0].rsplit("/", 1)[-1]
        self.copyright = copyright

        # set up module list
        self.modules = modules
        if module_env:
            module_envs = [m for m in os.environ.get(module_env, "").split(":") if m]
            if env_takes_precedence:
                self.modules = module_envs + self.modules
            else:
                self.modules += module_envs

    def autocomplete(self, tokens: list[str]):

        last_token = tokens[-1]

        # prepare the completion lists
        cmds = self.get_command_list()

        if len(tokens) == 1:
            completions = sorted(cmds)

        else:
            completions = sorted([opt for opt in cmds if opt.startswith(last_token)])

        # send results through fd 8
        ifs = os.environ.get("IFS", "\013")
        out_stream = os.fdopen(8, "w")
        out_stream.write(ifs.join(completions))
        sys.exit(0)

    def generic_usage(self):
        if self.copyright:
            _cerr(self.copyright)
        _cexit(
            f"  usage: {self.prog_name} [-i] [-l] [-h]\n"
            f"      {self.prog_name} CMD [ARGS]      run command or script (with ext .py)\n"
            f"      {self.prog_name} -i              run in interactive mode\n"
            f"      {self.prog_name} -l              list available commands\n"
            f"      {self.prog_name} -h              show this help message\n"
            f"  try: {self.prog_name} -l"
        )

    def generic_greet(self):
        _cerr(f"{self.prog_name} - command line interface")
        _cerr(f"Host: {platform.uname().node}")

    def generic_help(self):
        self.generic_usage()

    def version(self):
        _cerr(f"subcommand version: {__version__}")
        _cerr(f"try version command: {self.prog_name} version")

    def get_command_list(self) -> list[str]:
        # read sqpy.cmds directory
        cmd_files = []
        for module in self.modules:
            M = importlib.import_module(module)
            cmd_directory = pathlib.Path(M.__path__[0])
            cmd_files += cmd_directory.iterdir()
        cmds = []
        for p in cmd_files:
            if p.suffix != ".py":
                continue
            if p.name.startswith("_"):
                continue
            if "-" in p.name:
                _cerr(
                    "Warning: found command module filename containing dash: "
                    f"{p.name}"
                )
                continue
            cmds.append(p.stem.replace("_", "-"))
        return sorted(set(cmds))

    def show_commands(self):
        _cout("Available commands:")
        for cmd in self.get_command_list():
            _cout(f"  {cmd}")
        sys.exit(0)

    def run_main(
        self, main: Callable | None, init_argparser: Callable | None, args: list[str]
    ):

        if init_argparser is not None:
            parser = init_argparser()
            if not isinstance(parser, argparse.ArgumentParser):
                _cexit("ERR: init_argparser() did not return ArgumentParser instance")
            add_debug_to_parser(parser)

            # perform bash autocompletion if needed
            argcomplete.autocomplete(parser)

            # if autocomplete does not exit:
            args = parser.parse_args(args[1:])

            # provide a way to access the argument parser
            args.__parser__ = parser

            if main is not None:
                if args.debug:
                    from ipdb import launch_ipdb_on_exception

                    with launch_ipdb_on_exception():
                        _cerr("WARN: running in debug mode")
                        main(args)
                else:
                    main(args)
            else:
                _cexit("ERR: init_argparser() exists but main() does not")

        elif main is not None:
            import inspect

            if "args" in inspect.signature(main).parameters:
                main(args=args)
            else:
                main()

        # warn the user that nothing is being done
        else:
            raise RuntimeError(
                "FATAL ERROR: SubCommands.run_main() must be "
                "supplied with either main or init_argparser "
                "argument"
            )

    def run_cmd(self, args: list[str]):
        """this method will run module init_argparser() and main()"""

        command = args[0].replace("-", "_")
        L.debug("searching for module: %s", command)
        for module in self.modules:
            try:
                module_name = f"{module}.{command}"
                M = importlib.import_module(module_name)
                break
            except ModuleNotFoundError as err:
                if err.name != module_name:
                    raise
        else:
            _cexit(
                f"Command: {args[0]} does not exist! "
                f'Please check with "{self.prog_name} -l".'
            )

        _cerr(f"Executing: {M.__file__}")
        self.run_main(
            getattr(M, "main", None), getattr(M, "init_argparser", None), args
        )

    def run_script(self, args: list[str]):
        """this method will run an arbitrary python script"""

        path = args[0]
        _cerr(f"Executing: {path}")

        # expand home directory
        if path.startswith("~"):
            path = pathlib.Path(path).expanduser()
        else:
            path = pathlib.Path(path)

        if not path.exists():
            _cexit(f"Script: {path} does not exist!")

        with open(path) as fh:
            code = compile(fh.read(), path, "exec")
            exec(code, globals())
            _g = globals()
            _g["__name__"] = "__anyscript_main__"
            self.run_main(_g.get("main", None), _g.get("init_argparser", None), args)

    def main(self):

        tokens = []

        if "_ARGCOMPLETE" in os.environ:
            # we are invoked from argcomplete
            line = os.environ.get("COMP_LINE", "")
            tokens = line.split()
            if len(tokens) == 1 or (len(tokens) == 2 and not line.endswith(" ")):
                self.autocomplete(tokens)
            os.environ["COMP_LINE"] = line.split(" ", 1)[1]
            os.environ["COMP_POINT"] = str(len(os.environ["COMP_LINE"]))
            cmd = tokens[1]
        else:
            self.greet()
            if len(sys.argv) == 1:
                self.usage()
            cmd = sys.argv[1]

        if cmd == "-i" and self.allow_shell:
            # go to interactive mode
            import IPython

            IPython.embed()

        elif cmd == "-l":
            # show available commands
            self.show_commands()

        elif cmd == "-h":
            # show some help
            self.help()

        elif cmd == "-v":
            # show version
            self.version()

        elif cmd.endswith(".py") and self.allow_any_script:

            # run script, with:
            #   path = cmd
            #   args = sys.argv[2:] (eg. 'my_prog my_script arg0 arg1')
            self.run_script(tokens[1:] if any(tokens) else sys.argv[1:])

        else:
            self.run_cmd(tokens[1:] if any(tokens) else sys.argv[1:])


# EOF
