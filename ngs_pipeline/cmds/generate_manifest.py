__copyright__ = """
generate_manifest.py - ngs-pipeline command line
[https://github.com/vivaxgen/ngs-pipeline]

(c) 2022-2023 Hidayat Trimarsanto <trimarsanto@gmail.com>

All right reserved.
This software is licensed under MIT license.
Please read the README.txt of this software.
"""

# to improve the responsiveness during bash autocomplete, do not import heavy
# modules (such as numpy, pandas, etc) here, but instead import them within the
# functions that requires respective heavy modules

import sys
import os
import pathlib
from ngs_pipeline import cexit, cerr, arg_parser


def init_argparser():
    p = arg_parser("generate sample manifest file")
    p.add_argument("-o", "--outfile", default="outfile.tsv")

    m = p.add_mutually_exclusive_group()
    m.add_argument(
        "-s",
        "--single",
        default=False,
        action="store_true",
        help="fastq files are single (non-paired) such as ONT reads",
    )
    m.add_argument(
        "-p",
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

    p.add_argument(
        "--pause",
        type=int,
        default=0,
        help="pause (in seconds) before back to shell prompt, "
        "useful for automatic or batch processing so users can "
        "double check the sample names",
    )
    p.add_argument(
        "--ask-confirmation",
        default=False,
        action="store_true",
        help="ask confirmation to continue saving to output file",
    )

    p.add_argument(
        "-i",
        "--initial-manifest",
        default=None,
        type=str,
        help="optional initial manifest file to be appended",
    )

    p.add_argument("infiles", nargs="*")
    return p


def generate_manifest(args):

    # import heavy modules here if required
    import pandas as pd
    from ngs_pipeline import fileutils

    initial_df = None
    if args.initial_manifest:
        cerr(f"[Reading initial manifest file: {args.initial_manifest}]")
        initial_df = pd.read_table(args.initial_manifest)
        if not (("SAMPLE" in initial_df.columns) and ("FASTQ" in initial_df.columns)):
            cexit("[Initial manifest file does not have SAMPLE and/or FASTQ columns]")

    cerr(f"[Receiving {len(args.infiles)} source files]")
    if (initial_df is None) and (len(args.infiles) == 0):
        cexit("[No source files as input, please check your input files]")

    if args.single:
        mode = fileutils.ReadMode.SINGLETON
    elif args.paired:
        mode = fileutils.ReadMode.PAIRED_END
    else:
        mode = None

    sample_series = []
    fastq_series = []
    counter = 0

    if any(args.infiles):
        read_files = fileutils.ReadFileDict(
            args.infiles,
            args.underscore,
            args.remove_underscore_prefix,
            args.remove_prefix,
            mode,
        )

        if any(read_files.err_files):
            cexit(
                "ERROR: invalid input files: \n"
                + "\n".join(f"  {errmsg}" for errmsg in read_files.err_files)
            )

        # for each sample, process manifest line

        for sample in read_files.samples():
            sample_series.append(sample)
            items = [
                ",".join(item) if type(item) == tuple else item
                for item in read_files[sample]
            ]
            fastq_series.append(";".join(items))

    elif initial_df is None:
        cexit(
            "[No source files as input and no initial manifest file, please check both]"
        )

    df = pd.DataFrame(dict(SAMPLE=sample_series, FASTQ=fastq_series))
    if initial_df is not None:
        df = pd.concat([initial_df, df])
        # drop duplicate row
        df = df.drop_duplicates(keep="last")
        # check for duplicate sample name
        duplicate_samples = df.SAMPLE[df.SAMPLE.duplicated()]
        if any(duplicate_samples):
            cexit(
                "ERROR: duplicate sample code with differing read files"
                " found during merging with initial manifest:\n"
                + " ".join(list(duplicate_samples)),
                err_code=101,
            )

    # print at least 5 rows
    cerr("[Showing snippets of sample(s):]")
    cerr(str(df))
    cerr(
        "If the sample code is not correct, try to use option --underscore, "
        "--remove-underscore-prefix, or --remove-prefix"
    )

    if args.pause > 0:
        import time

        cerr(f"[Pausing for {args.pause} second(s) for manual inspection]")
        cerr("[Press CTRL-C to abort saving the output file]")
        time.sleep(args.pause)

    if args.ask_confirmation:
        resp = input(f"Continue saving to {args.outfile} [y/n]: ")
        if resp.lower().strip()[0] != "y":
            return

    df.to_csv(args.outfile, sep="\t", index=False)
    cerr(f"[Writing {len(df)} sample manifest to {args.outfile}]")


def main(args):
    generate_manifest(args)


# EOF
