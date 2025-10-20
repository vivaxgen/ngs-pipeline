__copyright__ = """
amtofastq.py - ngs-pipeline command line
[https://github.com/vivaxgen/ngs-pipeline]

(c) 2022-2024 Hidayat Trimarsanto <trimarsanto@gmail.com>

All right reserved.
This software is licensed under MIT license.
Please read the README.txt of this software.
"""

# to improve the responsiveness during bash autocomplete, do not import heavy
# modules (such as numpy, pandas, etc) here, but instead import them within the
# functions that requires respective heavy modules

import sys
import os

import ngs_pipeline
from ngs_pipeline import cerr, cexit, arg_parser, snakeutils


# usage: run_amtofastq.py


def init_argparser():
    p = snakeutils.init_argparser(
        "convert alignment map file (SAM/BAM/CRAM) to fastq files"
    )

    # p.add_argument('-j', type=int, default=72)

    p.add_argument(
        "--srcdir",
        default=".",
        help="Set source directory, where all the source files reside",
    )
    p.add_argument(
        "-s", default=None, help="Sample code, to be used for prefix of fastq filenames"
    )
    p.add_argument(
        "-i", default=None, help="Alignment map file to be splitted to fastq files"
    )
    p.add_argument(
        "infile",
        default=None,
        nargs="?",
        help="A tab-separated value with header SAMPLE and SOURCE",
    )
    return p


def amtofastq(args):

    from pathlib import Path

    # set SOURCE data

    if args.s:
        if not args.i:
            cexit("Please provide -i if you use -s")
        sources = {args.s: (Path(args.srcdir) / args.i).as_posix()}

    elif args.infile:
        import pandas as pd

        df = pd.read_table(args.infile)
        sources = {}
        for idx, row in df.iterrows():
            sources[row["SAMPLE"]] = (Path(args.srcdir) / row["SOURCE"]).as_posix()

    else:
        cexit("Please either provide INFILE or use -s & -i")

    # run smk

    config = dict(
        SOURCES=sources,
    )

    args.snakefile = snakeutils.get_snakefile_path(
        "amtofastq.smk", from_module=ngs_pipeline
    )

    status, elapsed_time = snakeutils.run_snakefile(args, config=config)

    if not status:
        cerr("[ERR: amtofastq did not successfully complete]")
    cerr(f"[Finish full run of amtofastq (time: {elapsed_time})]")


def main(args):
    amtofastq(args)


# EOF
