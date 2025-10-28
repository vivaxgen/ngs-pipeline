__copyright__ = """
run_targeted_variant_caller.py - ngs-pipeline command line
[https://github.com/vivaxgen/ngs-pipeline]

(c) 2023 Hidayat Trimarsanto <trimarsanto@gmail.com>

All right reserved.
This software is licensed under MIT license.
Please read the README.txt of this software.
"""

# to improve the responsiveness during bash autocomplete, do not import heavy
# modules (such as numpy, pandas, etc) here, but instead import them within the
# functions that requires respective heavy modules

import os
from ngs_pipeline import cerr, cexit, check_multiplexer
from ngs_pipeline.cmds import run_snakefile


def init_argparser():
    p = run_snakefile.init_argparser(desc="run targeted variant calling")
    p.arg_dict["snakefile"].choices = [
        "msf_targeted_varcall.smk",
        "panel_varcall_pe.smk",
        "panel_varcall_lr.smk",
        "msf_panel_varcall_pe.smk",
        "msf_panel_varcall_lr.smk",
    ]
    p.arg_dict["snakefile"].default = "msf_targeted_varcall.smk"

    # input/output options

    p.add_argument(
        "-o", "--outdir", default="analysis", help="directory for output [analysis/]"
    )

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
        "-u",
        "--underscore",
        type=int,
        default=0,
        help="number of consecutive underscore to be stripped from"
        "filenames to form sample code, counted in reverse",
    )

    p.add_argument(
        "--remove-underscore-prefix",
        type=int,
        default=0,
        help="the number of underscore to remove from the beginning of the filename",
    )

    p.add_argument(
        "--remove-prefix", default=None, help="prefix to remove from original filename"
    )

    p.add_argument("-i", "--manifest", default=None, help="manifest file  as input")
    p.add_argument(
        "infiles", nargs="*", help="FASTQ input files, eg. sample-1.fastq.gz"
    )

    return p


def run_targeted_variant_caller(args, optional_config={}):

    import pathlib
    import pickle
    from ngs_pipeline import fileutils

    # check exixtence of output directory
    if pathlib.Path(args.outdir).exists():
        while True:
            resp = input(
                f"Output directory {args.outdir} already exists.\n"
                f"Continue/overwrite/abort [c/o/a]: "
            )
            match resp.lower():

                case "o":
                    import shutil

                    shutil.rmtree(args.outdir)
                    break

                case "a":
                    cexit("Abort process.", err_code=101)

                case "c":
                    break

                case _:
                    cexit("Abort process.", err_code=101)

    # check we are inside a terminal multiplexer
    check_multiplexer(prompt=True)

    # get default snakefile from environment if none is provided
    if not args.snakefile:
        args.snakefile = os.environ.get("DEFAULT_SNAKEFILE", None)
        if args.snakefile:
            cerr(f"Obtaining snakefile from SNAKEFILE environment: {args.snakefile}")

    if not (any(args.infiles) or args.manifest):
        cexit(f"ERROR: need to have infiles or manifest file (--manifest)")

    if args.single:
        mode = fileutils.ReadMode.SINGLETON
    elif args.paired:
        mode = fileutils.ReadMode.PAIRED_END
    else:
        mode = None

    read_files = fileutils.ReadFileDict(
        args.infiles,
        underscore=args.underscore,
        underscore_prefix=args.remove_underscore_prefix,
        remove_prefix=args.remove_prefix,
        mode=mode,
        skip_list=[],
        manifest_file=args.manifest,
        sort_by_size=True,
    )

    # save the infiles to a pickle

    manifest_picklefile = pathlib.Path(args.outdir) / "metafile" / "manifest.pickle"
    manifest_picklefile.parent.mkdir(parents=True, exist_ok=True)

    with open(manifest_picklefile, "wb") as fout:
        pickle.dump(read_files, fout)
    cerr(f"Manifest file pickled for {len(read_files._d)} sample(s).")

    config = (
        dict(
            infiles=[],
            underscore=args.underscore,
            outdir=args.outdir,
            manifest_picklefile=manifest_picklefile,
        )
        | optional_config
    )
    status, elapsed_time = run_snakefile.run_snakefile(args, config=config)

    if not status:
        cerr("[WARNING: targeted variant calling did not successfully complete]")
    cerr(f"[Finish targeted variant calling (time: {elapsed_time})]")


def main(args):
    run_targeted_variant_caller(args)


# EOF
