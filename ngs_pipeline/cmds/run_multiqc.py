# run_multiqc.py - ngs-pipeline command
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2025, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

# to improve the responsiveness during bash autocomplete, do not import heavy
# modules (such as numpy, pandas, etc) here, but instead import them within the
# functions that require the respective heavy modules

import os
import pathlib
from ngs_pipeline import cerr, cexit, get_snakefile_path
from ngs_pipeline.cmds import run_snakefile

# this is a wrapper to run initialize.smk


def init_argparser():

    p = run_snakefile.init_argparser(
        "prepare and initialize all required files and settings"
    )
    p.add_argument(
        "--outdir",
        default="00-QC",
        help="output directory (relative to input directory unless absolute path is given)",
    )
    p.add_argument("indir", help="input directory containing FASTQ files")
    return p


def run_multiqc(args):

    import ngs_pipeline

    os.environ["NGSENV_BASEDIR"] = os.environ["NGS_PIPELINE_BASE"]

    args.snakefile = get_snakefile_path("multiqc.smk", from_module=ngs_pipeline)
    args.no_config_cascade = True
    args.force = True

    run_snakefile.run_snakefile(args, dict(indir=args.indir, outdir=args.outdir))


def main(args):
    run_multiqc(args)


# EOF
