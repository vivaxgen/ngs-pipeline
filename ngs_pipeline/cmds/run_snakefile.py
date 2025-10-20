#!/usr/bin/env ngs-pl
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(c) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"


import os
import sys
import pathlib
import ngs_pipeline
from ngs_pipeline import (
    cerr,
    check_NGSENV_BASEDIR,
    check_NGS_PIPELINE_BASE,
    check_force,
    setup_config,
    snakeutils,
    subcommands,
)


get_snakefile_path = snakeutils.get_snakefile_path

subcommands.set_argument_parser_class(snakeutils.ArgumentParser)


def init_argparser(desc=None):

    p = snakeutils.init_argparser(
        p=subcommands.arg_parser(desc=desc or "run arbitrary snakefile")
    )
    return p


def run_snakefile(
    args,
    config: dict = {},
    additional_cli_args="",
    workdir: str | pathlib.Path | None = None,
    show_status: bool = True,
    show_config_files: bool = False,
):

    env_base_dir = pathlib.Path(check_NGSENV_BASEDIR())
    pipeline_base_dir = pathlib.Path(check_NGS_PIPELINE_BASE())

    executor = snakeutils.SnakeExecutor(
        args,
        setup_config_func=setup_config,
        workdir=workdir,
        env_basedir=env_base_dir,
        from_module=ngs_pipeline,
        default_config_file=pipeline_base_dir / "etc" / "default-config.yaml",
        show_config_files=show_config_files,
    )

    status, elapsed_time = executor.run(
        snakefile=args.snakefile,
        config=config,
        additional_cli_args=additional_cli_args,
        force=check_force(args.force),
        no_config_cascade=snakeutils.check_env("NGS_PIPELINE_NO_CONFIG_CASCADE")
        or args.no_config_cascade,
    )

    if show_status:
        if not status:
            cerr(
                f"[WARNING: snakefile {args.snakefile} did not successfully "
                "complete]"
            )
        cerr(f"[Finish running snakefile {args.snakefile} (time: {elapsed_time})]")

    return status, elapsed_time


def main(args):
    run_snakefile(args)


# EOF
