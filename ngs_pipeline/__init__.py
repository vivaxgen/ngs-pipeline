# SPDX-FileCopyrightText: 2023-2006 Hidayat Trimarsanto <trimarsanto@gmail.com>
# SPDX-License-Identifier: MIT

from __future__ import annotations

__copyright__ = "(c) 2023-2006 Hidayat (Anto) Trimarsanto <trimarsanto@gmail.com>"
__author__ = "trimarsanto@gmail.com"
__license__ = "MIT"

# __init__.py - ngs-pipeline package
# [https://github.com/vivaxgen/ngs-pipeline]

# library providing file common functions

import os
import sys
import argparse
import platform
import argcomplete
import pathlib
import types

from .subcommands import arg_parser


def cout(msg):
    print(msg, file=sys.stdout)


def cerr(msg):
    print(msg, file=sys.stderr)
    sys.stderr.flush()


def cexit(msg, err_code=1):
    cerr(msg)
    sys.exit(err_code)


def greet():
    cerr(
        f'{sys.argv[0].split("/")[-1]} - ngs-pipeline command line interface\n'
        f"[https://github.com/vivaxgen/ngs-pipeline]"
    )
    cerr(f"Host: {platform.uname().node}")


def check_VVG_BASEDIR() -> str:
    if "VVG_BASEDIR" not in os.environ:
        cexit(
            "ERROR: VVG_BASEDIR environment is not set. "
            "Please set proper shell enviroment by sourcing relevant activate.sh"
        )
    return os.environ["VVG_BASEDIR"]


def check_NGSENV_BASEDIR() -> str:
    if "NGSENV_BASEDIR" not in os.environ:
        cexit(
            "ERROR: NGSENV_BASEDIR environment is not set. "
            "Please set proper shell enviroment by sourcing relevant activate.sh"
        )
    return os.environ["NGSENV_BASEDIR"]


def check_NGS_PIPELINE_BASE() -> str:
    if "NGS_PIPELINE_BASE" not in os.environ:
        cexit(
            "ERROR: NGS_PIPELINE_BASE environment is not set. "
            "Please set proper shell enviroment by sourcing relevant activate.sh"
        )
    return os.environ["NGS_PIPELINE_BASE"]


def check_multiplexer(prompt=False, remote_only=True) -> bool:
    if remote_only and not os.environ.get("SSH_TTY", None):
        return True
    if os.environ.get("NGS_IGNORE_TERM_MULTIPLEXER_CHECK", None):
        return True
    if "TMUX" in os.environ:
        return True
    if "screen" in os.environ.get("TERM", "") and "screen" in os.environ.get(
        "TERMCAP", ""
    ):
        return True
    if prompt:
        from ngs_pipeline.timed_input import timed_input

        resp = timed_input(
            "Warning! You are not in a terminal multiplexer interactive session.\n"
            "Still continue [Y/n]: ",
            10,
            "y",
        )
        if not resp.lower().startswith("y"):
            cerr(f"[Process terminated.]")
            sys.exit(101)
    return False


def check_force(force: bool) -> bool:
    return force or bool(os.environ.get("NGS_PIPELINE_FORCE", 0))


def add_pgline(alignment_file, pg_dict) -> dict:
    """append new PG line using a PG dict, return the full header"""

    header = alignment_file.header.to_dict()
    pg_header = header["PG"]
    pg_dict["PP"] = pg_header[-1]["ID"]
    pg_header.append(pg_dict)
    return header


def get_mode(filename, mode) -> str:
    """return either (r, rb, w, wb) depending on file name"""
    if filename.endswith(".bam"):
        return mode + "b"
    return mode


def setup_config(d={}) -> dict:
    d["NGSENV_BASEDIR"] = check_NGSENV_BASEDIR()
    d["NGS_PIPELINE_BASE"] = check_NGS_PIPELINE_BASE()

    def _path_to_str(value):
        if isinstance(value, pathlib.Path):
            return value.as_posix()
        if isinstance(value, dict):
            return {k: _path_to_str(v) for k, v in value.items()}
        if isinstance(value, list):
            return [_path_to_str(v) for v in value]
        if isinstance(value, tuple):
            return tuple(_path_to_str(v) for v in value)
        if isinstance(value, set):
            return {_path_to_str(v) for v in value}
        return value

    # this function needs to return a dict

    for k, v in d.items():
        d[k] = _path_to_str(v)

    return d


def prepare_command_log():
    import time

    return {
        "lastInvoked": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
        "command": " ".join(sys.argv),
        "NGS_PIPELINE_BASE": check_NGS_PIPELINE_BASE(),
        "NGSENV_BASEDIR": check_NGSENV_BASEDIR(),
    }


# EOF
