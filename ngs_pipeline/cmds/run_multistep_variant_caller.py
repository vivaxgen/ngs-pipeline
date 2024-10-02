# run_multistep_variant_caller.py - ngs-pipeline command
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

# to improve the responsiveness during bash autocomplete, do not import heavy
# modules (such as numpy, pandas, etc) here, but instead import them within the
# functions that require the respective heavy modules

import os
import pathlib
from ngs_pipeline import (
    cerr,
    cexit,
    arg_parser,
    get_snakefile_path,
    check_multiplexer,
    snakeutils,
)

# this is a wrapper to run sample_preparation, individual sample genotyping,
# and joint variant calling in a single command


def init_argparser():
    p = snakeutils.init_argparser("run whole steps of discovery variant caller")

    # file input arguments, similar to generate-manifest command
    m = p.add_mutually_exclusive_group()
    m.add_argument(
        "--single",
        default=False,
        action="store_true",
        help="fastq files are single (non-paired) such as ONT reads",
    )
    m.add_argument(
        "--paired",
        default=False,
        action="store_true",
        help="fastq files are paired such as Illumina paired-end reads",
    )
    p.add_argument(
        "-P",
        default=None,
        help="procfile containing number of task/sample to process in parallel",
    )

    p.add_argument(
        "-u",
        "--underscore",
        type=int,
        default=0,
        help="no of consecutive underscore to be stripped from "
        "filenames to form sample code, counted in reverse",
    )

    p.add_argument("-i", "--manifest", default=None, help="manifest file  as input")

    p.add_argument("-o", "--outdir", required=True, help="output directory")
    p.add_argument("infiles", nargs="*", help="read files")

    return p


def run_multistep_variant_caller(args, console=True):

    import ngs_pipeline
    from ngs_pipeline import get_snakefile_path, snakeutils
    import sys

    # check we are inside a terminal multiplexer
    check_multiplexer(prompt=True)

    # prevent other commands to check for terminal multiplexer
    os.environ["NGS_IGNORE_TERM_MULTIPLEXER_CHECK"] = "1"

    config = dict(
        manifest=args.manifest,
        underscore=args.underscore,
        singleton=args.single,
        paired_end=args.paired,
        outdir=args.outdir,
        infiles=args.infiles,
        jobs=args.j,
        procfile=args.P,
        rerun=(args.manifest is None) and not any(args.infiles),
        unlock=args.unlock,
    )

    args.snakefile = get_snakefile_path(
        "multistep_variant_calling.smk", from_module=ngs_pipeline
    )

    if config["rerun"]:
        cerr("[INFO: unlocking all working directories...]")
        # unlock for main snakefile and all necessary rules
        args.unlock = config["unlock"] = True
        # args.rerun = config['rerun'] = False

        # unlock for main snakefile
        status, elapse_time = snakeutils.run_snakefile(args, config=config)

        #  unlock for all necessary rules
        args.unlock = False
        status, elapse_time = snakeutils.run_snakefile(args, config=config)

        # no more unlock
        config["unlock"] = False
        # args.rerun = config['rerun'] = True
        cerr("[INFO: finished unlocking all working directories...]")

    status, elapsed_time = snakeutils.run_snakefile(args, config=config)
    if not console:
        return status, elapsed_time

    if not status:
        cerr(
            "[WARNING: full run of multistep variant calling "
            "did not successfully complete]"
        )
    cerr(f"[Finish full run of multi-step variant calling (time: {elapsed_time})]")


def main(args):
    run_multistep_variant_caller(args)


# EOF
