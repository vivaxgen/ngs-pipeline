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
    setup_config,
    snakeutils,
)


def init_argparser(desc=None):
    p = snakeutils.init_argparser(desc=desc or "run arbitrary snakefile")
    return p


def run_snakefile(args, config: dict = {}, show_status: bool = True):

    env_base_dir = pathlib.Path(check_NGSENV_BASEDIR())
    pipeline_base_dir = pathlib.Path(check_NGS_PIPELINE_BASE())

    executor = snakeutils.SnakeExecutor(
        args,
        setup_config_func=setup_config,
        env_basedir=env_base_dir,
        from_module=ngs_pipeline,
        default_config_file=pipeline_base_dir / "etc" / "default-config.yaml",
    )

    status, elapsed_time = executor.run(
        snakefile=args.snakefile,
        config=config,
        force=snakeutils.check_env("NGS_PIPELINE_FORCE") or args.force,
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
