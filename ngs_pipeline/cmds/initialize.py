# initialize.py - ngs-pipeline command
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

# to improve the responsiveness during bash autocomplete, do not import heavy
# modules (such as numpy, pandas, etc) here, but instead import them within the
# functions that require the respective heavy modules

import os
import pathlib
from ngs_pipeline import cerr, cexit
from ngs_pipeline.cmds import run_snakefile

# this is a wrapper to run initialize.smk


def init_argparser():

    p = run_snakefile.init_argparser(
        "prepare and initialize all required files and settings"
    )
    p.add_argument("--mem_gb", type=int, default=16, help="required memory (in MB)")
    return p


def initialize(args):

    import ngs_pipeline

    args.snakefile = run_snakefile.get_snakefile_path(
        "initialize.smk", from_module=ngs_pipeline
    )
    args.no_config_cascade = True
    args.force = True

    run_snakefile.run_snakefile(args, dict(index_mem_gb=args.mem_gb))


def main(args):
    initialize(args)


# EOF
