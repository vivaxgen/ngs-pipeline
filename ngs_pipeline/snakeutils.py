# SPDX-FileCopyrightText: 2024-2006 Hidayat Trimarsanto <trimarsanto@gmail.com>
# SPDX-License-Identifier: MIT

from __future__ import annotations

__copyright__ = "(c) 2024-2006 Hidayat (Anto) Trimarsanto <trimarsanto@gmail.com>"
__author__ = "trimarsanto@gmail.com"
__license__ = "MIT"

# snakeutils.py
# [https://github.com/trmznt/py-snakeutils]

__version__ = "2026.07.19.01"

# this module provides wrapper to execute Snakemake file from Python code


import os
import sys
import pathlib
import time
import datetime
import logging
import argparse
import types
import importlib

from abc import ABC, abstractmethod
from typing import Callable, Iterable, Iterator, Sequence, Any

L = logging.getLogger(__name__)


_DEFAULT_RULE_PATH = None


def _cout(msg: str) -> None:
    print(msg, file=sys.stdout)


def _cerr(msg: str) -> None:
    print(msg, file=sys.stderr)
    sys.stderr.flush()


def _cexit(msg: str, exit_code: int = 1) -> None:
    _cerr(msg)
    sys.exit(exit_code)


def set_default_rule_path(module: types.ModuleType, overwrite: bool = False) -> None:
    global _DEFAULT_RULE_PATH

    if _DEFAULT_RULE_PATH is not None and not overwrite:
        return

    if not hasattr(module, "__path__") or not module.__path__:
        raise ValueError(f"Module '{module.__name__}' does not have a valid __path__.")

    _DEFAULT_RULE_PATH = pathlib.Path(module.__path__[0]) / "rules"


#
# argument parser helpers
#

FilePredicate = Callable[[pathlib.Path, argparse.Namespace], bool]


class ValueProvider(ABC):
    """Provides values for completion or validation."""

    @abstractmethod
    def values(self) -> Sequence[str]:
        """Return available values."""
        raise NotImplementedError


class LazyProvider(ValueProvider):
    """
    Lazily computes values once and caches them.

    The cache can be invalidated explicitly.
    """

    def __init__(
        self,
        loader: Callable[[], Iterable[str]],
    ) -> None:
        self._loader = loader
        self._cache: list[str] | None = None

    def values(self) -> list[str]:
        if self._cache is None:
            self._cache = list(self._loader())
        return self._cache

    def invalidate(self) -> None:
        """Clear the cached values."""
        self._cache = None

    def reload(self) -> list[str]:
        """Reload and return the latest values."""
        self.invalidate()
        return self.values()


class StaticProvider(ValueProvider):
    """Provider backed by static values."""

    def __init__(
        self,
        values: Iterable[str],
    ) -> None:
        self._values = list(values)

    def values(self) -> Sequence[str]:
        return self._values


class CompletionSource(ABC):
    """Base completion source."""

    @abstractmethod
    def help_items(self) -> Sequence[str]: ...

    @abstractmethod
    def complete(
        self,
        prefix: str,
        parsed_args: object,
        **kwargs,
    ) -> Iterable[str]: ...


class SuggestionsSource(CompletionSource):
    """Completion from a ValueProvider."""

    def __init__(
        self,
        provider: ValueProvider,
        *,
        label: str = "Suggested:",
        help_limit: int = 8,
    ) -> None:

        self.provider = provider
        self.label = label
        self.help_limit = help_limit

    def help_items(self) -> Sequence[str]:

        values = list(self.provider.values())

        if len(values) > self.help_limit:
            return [
                *values[: self.help_limit],
                "...",
            ]

        return values

    def complete(
        self,
        prefix: str,
        parsed_args: object,
        **kwargs,
    ) -> Iterator[str]:
        for item in self.provider.values():
            if item.startswith(prefix):
                yield item


class FilesSource(CompletionSource):
    """
    Filesystem completion source.

    Parameters
    ----------
    predicate
        Optional callback used to decide whether a filesystem entry
        should be included.

        The callback receives

        - the candidate Path
        - the parsed argparse Namespace
    """

    def __init__(
        self,
        *,
        predicate: FilePredicate | None = None,
    ) -> None:

        self._predicate = predicate

    def help_items(self) -> list[str]:
        return []

    def complete(
        self,
        prefix: str,
        parsed_args: argparse.Namespace,
        **kwargs,
    ) -> Iterator[str]:

        # Preserve whether user typed:
        #   workflows/
        # versus
        #   workflows
        has_trailing_slash = prefix.endswith("/")

        if has_trailing_slash:

            parent = pathlib.Path(prefix)
            partial = ""

        else:

            path = pathlib.Path(prefix)
            parent = (
                path.parent if path.parent != pathlib.Path("") else pathlib.Path(".")
            )
            partial = path.name

        try:
            children = sorted(
                parent.iterdir(),
                key=lambda p: p.name,
            )

        except OSError:
            return

        for child in children:

            if partial and not child.name.startswith(partial):
                continue

            if self._predicate is not None:

                try:
                    if not self._predicate(
                        child,
                        parsed_args,
                    ):
                        continue

                except Exception:
                    continue

            if child.is_dir():
                yield str(child) + "/"
            else:
                yield str(child)


class Argument:
    """
    Wrapper around argparse.Action.

    Unknown attributes are delegated to argparse.Action.
    """

    def __init__(
        self,
        action: Any,
    ) -> None:

        object.__setattr__(
            self,
            "_action",
            action,
        )

        object.__setattr__(
            self,
            "_completion_sources",
            [],
        )

        action._argument_wrapper = self

    def __getattr__(
        self,
        name: str,
    ) -> Any:

        return getattr(
            self._action,
            name,
        )

    def __setattr__(
        self,
        name: str,
        value: Any,
    ) -> None:

        if name.startswith("_"):
            object.__setattr__(
                self,
                name,
                value,
            )
        else:
            setattr(
                self._action,
                name,
                value,
            )

    @property
    def completion_sources(
        self,
    ) -> list[CompletionSource]:

        return self._completion_sources

    def suggestions(
        self,
        provider: ValueProvider | Callable[[], Iterable[str]],
        *,
        label: str = "Suggested",
        help_limit: int = 8,
    ) -> Argument:
        """
        Configure suggestion completion.

        Parameters
        ----------
        provider:
            ValueProvider or callable returning strings.

        label:
            Text shown in --help.

        help_limit:
            Number of examples displayed.
        """

        if not isinstance(
            provider,
            ValueProvider,
        ):
            provider = LazyProvider(provider)

        self._completion_sources[:] = [
            source
            for source in self._completion_sources
            if not isinstance(
                source,
                SuggestionsSource,
            )
        ]

        self._completion_sources.append(
            SuggestionsSource(
                provider,
                label=label,
                help_limit=help_limit,
            )
        )

        return self

    def files(
        self,
        *,
        predicate: FilePredicate | None = None,
    ) -> Argument:
        """
        Enable filesystem completion.

        Parameters
        ----------
        predicate
            Optional callback that filters filesystem entries.

            example:

            predicate=lambda path, args: (
                path.is_dir()
                or path.suffix == ".vcf"
            )

        """

        self._completion_sources[:] = [
            s for s in self._completion_sources if not isinstance(s, FilesSource)
        ]

        self._completion_sources.append(
            FilesSource(
                predicate=predicate,
            )
        )

        return self

    def complete(
        self,
        prefix: str,
        parsed_args: object,
        **kwargs: Any,
    ) -> Iterator[str]:

        seen = set()

        for source in self.completion_sources:
            try:
                for item in source.complete(
                    prefix,
                    parsed_args,
                    **kwargs,
                ):
                    if item not in seen:
                        seen.add(item)
                        yield item
            except Exception:
                pass


class SuggestionHelpFormatter(argparse.HelpFormatter):
    """
    Help formatter that displays completion previews.
    """

    def _expand_help(
        self,
        action: argparse.Action,
    ) -> str:

        text = super()._expand_help(action)

        argument = getattr(
            action,
            "_argument_wrapper",
            None,
        )

        if argument is None:
            return text

        for source in argument.completion_sources:

            try:
                items = source.help_items()

            except Exception:
                continue

            if not items:
                continue

            label = getattr(
                source,
                "label",
                "Available",
            )

            text += (
                "\n" + " " * self._current_indent + f"{label}: " + "\n, ".join(items)
            )

        return text


class ArgumentParser(argparse.ArgumentParser):
    """
    Extended argparse parser.

    Features
    --------
    - Keeps argparse behavior unchanged.
    - Provides register_argument().
    - Maintains arg_dict externally.
    - Supports argcomplete integration.
    """

    LazyProvider = LazyProvider

    def __init__(
        self,
        *args: Any,
        **kwargs: Any,
    ) -> None:

        kwargs.setdefault(
            "formatter_class",
            SuggestionHelpFormatter,
        )

        super().__init__(
            *args,
            **kwargs,
        )

        self.arg_dict: dict[str, Argument] = {}

    def register_argument(
        self,
        *args: Any,
        **kwargs: Any,
    ) -> Argument:
        """
        Add an argument and return an Argument wrapper.
        """

        action = super().add_argument(
            *args,
            **kwargs,
        )

        return Argument(action)

    def enable_completion(self) -> None:
        """
        Enable argcomplete support.

        This must be called after all arguments are registered.
        """

        import argcomplete

        for action in self._actions:

            argument = getattr(
                action,
                "_argument_wrapper",
                None,
            )

            if argument is None:
                continue

            def completer(
                prefix: str,
                parsed_args: object,
                argument: Argument = argument,
                **kwargs: Any,
            ):
                return list(
                    argument.complete(
                        prefix,
                        parsed_args,
                        **kwargs,
                    )
                )

            action.completer = completer  # type: ignore

        argcomplete.autocomplete(self)


__ARGUMENT_PARSER__ = ArgumentParser


def set_argument_parser_class(class_) -> None:
    global __ARGUMENT_PARSER__
    __ARGUMENT_PARSER__ = class_


def init_argparser(desc: str = "", p: ArgumentParser | None = None) -> ArgumentParser:
    """provide common arguments for snakemake-based cli"""

    p = p or __ARGUMENT_PARSER__(description=desc)
    # p.arg_dict = {}

    # snakemake arguments
    p.add_argument(
        "-j", type=int, default=32, help="number of jobs to be executed in parallel"
    )

    # debugging process
    p.add_argument(
        "--dry-run", default=False, action="store_true", help="run in dry-mode"
    )
    p.add_argument(
        "--showcmds",
        default=False,
        action="store_true",
        help="show shell commands, similar to --printshellcmds",
    )
    p.add_argument(
        "-p",
        "--printshellcmds",
        default=False,
        action="store_true",
        help="show shell commands",
    )
    p.add_argument(
        "--show-config-files",
        default=False,
        action="store_true",
        help="show the order of config files to parse",
    )
    p.add_argument(
        "--keep-incomplete",
        default=False,
        action="store_true",
        help="keep incomplete files",
    )

    # continuation of previous run
    p.add_argument(
        "--unlock",
        default=False,
        action="store_true",
        help="unlock the directory from unfinished/failed snakemake run",
    )
    p.add_argument(
        "--rerun",
        default=False,
        action="store_true",
        help="continue running the snakemake workflow from previoius point",
    )
    p.add_argument(
        "--touch",
        default=False,
        action="store_true",
        help="touch all output files, to avoid re-running the ngs-pipeline "
        "such as after modifying/debugging snakemake file",
    )

    # running configuration
    p.add_argument(
        "--profile",
        default=None,
        help="snakemake profile to be used, also can be set from SNAKEMAKE_PROFILE env",
    )
    p.add_argument(
        "--nocluster",
        default=False,
        action="store_true",
        help="run without cluster support (eg. only on local node), useful for debugging",
    )

    # general options
    p.arg_dict["target"] = p.register_argument(
        "-t",
        "--target",
        default=[],
        action="append",
        help="target rule(s) in the snakefile [all]",
    )
    p.arg_dict["snakefile"] = p.register_argument(
        "--snakefile", default=None, help="snakemake file to be called"
    )

    # configuration files
    p.add_argument(
        "--base-config",
        default=None,
        help="path for base configuration file, relative to "
        "base environment directory",
    )
    p.arg_dict["panel"] = p.register_argument(
        "--panel",
        default=None,
        help="panel to be used (eg. PANEL -> configs/PANEL.yaml as base config)",
    )
    p.add_argument(
        "-c", "--config", default=[], action="append", help="config file(s) to append"
    )
    p.add_argument(
        "-f",
        "--force",
        default=False,
        action="store_true",
        help="force the processing even if the working directory is not "
        "under current pipeline environment base directory",
    )
    p.add_argument(
        "--no-config-cascade",
        default=False,
        action="store_true",
        help="prevent from reading cascading configuration file",
    )

    return p


def check_env(env_name: str) -> bool:
    if env_name in os.environ:
        return True
    return False


def setup_config(config) -> dict:
    # dummy setup configurator
    return config


def path_to_str(value) -> Any:
    if isinstance(value, pathlib.Path):
        return value.as_posix()
    if isinstance(value, dict):
        return {k: path_to_str(v) for k, v in value.items()}
    if isinstance(value, list):
        return [path_to_str(v) for v in value]
    if isinstance(value, tuple):
        return tuple(path_to_str(v) for v in value)
    if isinstance(value, set):
        return {path_to_str(v) for v in value}
    return value


class SnakeExecutor(object):

    def __init__(
        self,
        # arguments parsed from the init_argparser above
        args,
        *,
        # basic configuration
        setup_config_func: Callable[
            [
                dict,
            ],
            dict,
        ] = setup_config,
        # working directory
        workdir: str | pathlib.Path | None = None,
        # show configuration files in the parsing order to stderr
        show_config_files: bool = False,
        # environment base dir as root for cascading configuration
        env_basedir: str | pathlib.Path | None = None,
        # module to get the snakemake file from
        from_module: types.ModuleType | None = None,
        # default configuration file (eg from software installation)
        # if exists, this would be added as the first configuration file to parse
        default_config_file: str | pathlib.Path | None = None,
    ) -> None:

        from snakemake import cli

        self.args = args
        self.setup_config_func = setup_config_func
        self.workdir = pathlib.Path(workdir) if workdir else None
        self.show_config_files = show_config_files or args.show_config_files
        self.env_basedir = (
            pathlib.Path(env_basedir) if env_basedir else pathlib.Path.cwd()
        )
        self.from_module = from_module
        self.default_config_file = (
            pathlib.Path(default_config_file) if default_config_file else None
        )

        if self.from_module:
            set_default_rule_path(self.from_module)

        # monkey-patch snakemake cli class
        cli_parse_config = cli.parse_config

        def parse_config(entries):
            if type(entries) == dict:
                return entries
            return cli_parse_config(entries)

        cli.parse_config = parse_config
        # end of monkey patching

    def run(
        self,
        # snakefile to run
        snakefile: str | pathlib.Path | None = None,
        # configuration to append/update
        config: dict = {},
        additional_cli_args: str = "",
        *,
        from_module: types.ModuleType | None = None,
        # allow to run the snakefile even if not inside environment directory
        force: bool = False,
        # prevent from performing cascading configuration parsing
        no_config_cascade: bool = False,
    ) -> tuple[int, datetime.timedelta]:

        from snakemake import cli
        import shlex

        cwd = self.workdir or pathlib.Path.cwd()
        if "__workdir__" in config:
            raise ValueError('ERR: config key "__workdir__" is reserved')
        config["__workdir__"] = cwd
        _cerr(f"Current working directory: {cwd}")

        if not (force or self.args.force) and not cwd.is_relative_to(self.env_basedir):
            _cexit(
                f"ERROR: current directory {cwd} is not relative to {self.env_basedir}"
            )

        snakefile = get_snakefile_path(snakefile or self.args.snakefile)

        # check sanity
        if not snakefile:
            _cexit(
                "ERR: Please provide snakefile to execute using --snakefile argument."
            )

        # process config files and add to config __configfiles__ key

        configfiles = [pathlib.Path(cf) for cf in reversed(self.args.config)]
        if "__configfiles__" in config:
            raise ValueError('ERR: config key "__configfiles__" is reserved')
        config["__configfiles__"] = configfiles

        if no_config_cascade or self.args.no_config_cascade:
            if (configfile := self.env_basedir / "config.yaml").is_file():
                configfiles.append(configfile)
            config_dirs = [cwd]
        else:
            L.debug("processing cascading configuration")
            # for each config directory, check config file existence
            config_dirs = []
            config_path = cwd
            while config_path.is_relative_to(self.env_basedir):
                config_dirs.append(config_path)
                configfile = config_path / "config.yaml"
                if configfile.is_file():
                    configfiles.append(configfile)
                config_path = config_path.parent

        # get panel configuration and set as base configuration from configs/
        if self.args.panel:
            if self.args.base_config:
                _cexit(f"ERROR: cannot use both --panel and --base-config")
            self.args.base_config = "configs/" + self.args.panel + ".yaml"

        if self.args.base_config:
            if is_abs_or_rel_path(self.args.base_config):
                # real full path
                configfiles.append(self.args.base_config)
            else:
                configfiles.append(self.env_basedir / self.args.base_config)

        # get config file at root of environment base directory
        # this config file can be overridden by base/panel/custom config files
        if (root_configfile := self.env_basedir / "configs" / "config.yaml").is_file():
            configfiles.append(root_configfile)

        # if provided, this will be the default config file that comes with
        # the software/tool installation
        if self.default_config_file:
            configfiles.append(self.default_config_file)

        if not any(configfiles):
            _cexit(f"ERROR: cannot find any config.yaml in {config_dirs}")

        configfiles.reverse()
        if self.show_config_files:
            _cerr(f"config_files: {configfiles}")

        # setting up profile

        if "SNAKEMAKE_PROFILE" in os.environ:
            if self.args.profile is None:
                self.args.profile = os.environ["SNAKEMAKE_PROFILE"]

        if self.args.nocluster or "SNAKEMAKE_NOCLUSTER" in os.environ:
            # nocluster means prevent from running using batch/job scheduler,
            # then use special profile "none"
            self.args.profile = "none"

        # set targets
        if type(self.args.target) != list:
            targets = [self.args.target]
        else:
            targets = self.args.target if any(self.args.target) else ["all"]

        try:
            # set with --profile first
            if self.args.profile:
                argv = ["--profile", self.args.profile]
            else:
                argv = []

            # add compatible additional arguments from command line
            if self.args.keep_incomplete:
                argv.append("--keep-incomplete")

            # extend the arguments with additional arguments from the function call
            argv.extend(shlex.split(additional_cli_args))

            # extra argument pass from SNAKEMAKE_EXTRA_ARGS environment variable
            if "SNAKEMAKE_EXTRA_ARGS" in os.environ:
                argv.extend(shlex.split(os.environ["SNAKEMAKE_EXTRA_ARGS"]))

            # XXX: need to modify to use snakemake API
            L.debug("parsing snakemake arguments")
            parser, args = cli.parse_args(argv)

            args.snakefile = get_snakefile_path(
                snakefile, from_module=from_module or self.from_module
            )

            # set args further from self.args
            args.configfile = configfiles
            # cargs.config = [f'{k}={v}' for k, v in setup_config(config).items()]
            args.config = path_to_str(self.setup_config_func(config))
            args.targets = targets

            # running mode
            args.dryrun = self.args.dry_run
            args.reason = True if self.args.dry_run else False
            args.keep_incomplete = self.args.keep_incomplete
            args.touch = self.args.touch
            args.rerun_incomplete = self.args.rerun
            args.unlock = self.args.unlock

            # running parameters
            args.cores = self.args.j
            args.printshellcmds = self.args.showcmds or self.args.printshellcmds

            L.debug("invoking snakemake client")
            _cerr(f"Running snakefile: {args.snakefile}")
            start_time = time.monotonic()
            status = cli.args_to_api(args, parser)
            finish_time = time.monotonic()

            return (status, datetime.timedelta(seconds=int(finish_time - start_time)))

        except Exception as e:
            cli.print_exception(e)
            sys.exit(1)

        raise RuntimeError("FATAL ERROR: should not execute this part of code")


def get_snakefile_path(
    filepath: str | pathlib.Path,
    snakefile_root: pathlib.Path | None = None,
    from_module: types.ModuleType | None = None,
    strict_mode: bool = True,
) -> pathlib.Path:
    """
    - return real path of snakefile
    - filepath can be string or pathlib.Path with either absolute, relative or plain filename
    - filepath can also be in the format of module::filename where it is expected to have
      module/rules/filename path structure.
    """

    if type(filepath) == str and "::" in filepath:
        module_name, filepath = filepath.split("::")
        from_module = importlib.import_module(module_name)
    elif is_abs_or_rel_path(filepath):
        if type(filepath) == str:
            return pathlib.Path(filepath)
        return filepath  # type: ignore

    if from_module is not None:
        snakefile_root = pathlib.Path(from_module.__path__[0]) / "rules"
    if snakefile_root is not None:
        return snakefile_root / filepath
    if _DEFAULT_RULE_PATH and strict_mode:
        return _DEFAULT_RULE_PATH / filepath
    if not strict_mode:
        return pathlib.Path(filepath)
    raise ValueError(f"ERR: cannot determine the path for {filepath}")


def is_abs_or_rel_path(filepath: str | pathlib.Path) -> bool:
    filepath = filepath.as_posix() if isinstance(filepath, pathlib.Path) else filepath
    if (
        filepath.startswith("/")
        or filepath.startswith("./")
        or filepath.startswith("../")
    ):
        return True
    return False


def scan_for_config_keywords(path) -> list[str]:
    """return a list of keywords used as config keys in any of the
    snakemake rules
    """

    import re

    mo_bracket = re.compile(r"""config\s*\[\s*['"]([^]]*)['"]\s*\]""")
    mo_get = re.compile(r"""config\s*\.\s*get\s*\(\s*['"]([^'"]*)""")

    keywords = []

    path = pathlib.Path(path)

    for rule_file in path.glob("*.smk"):
        with open(rule_file) as f_in:
            for line in f_in:
                keys = mo_bracket.findall(line) + mo_get.findall(line)
                if any(keys):
                    keywords += keys

    return sorted(set(keywords))


# EOF
