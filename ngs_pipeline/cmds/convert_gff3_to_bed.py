# consolidate_samples.py - ngs-pipeline command
# [https://github.com/vivaxgen/ngs-pipeline]

__author__ = "Hidayat Trimarsanto"
__copyright__ = "(C) 2025, Hidayat Trimarsanto"
__email__ = "trimarsanto@gmail.com,hidayat.trimarsanto@menzies.edu.au"
__license__ = "MIT"

# to improve the responsiveness during bash autocomplete, do not import heavy
# modules (such as numpy, pandas, etc) here, but instead import them within the
# functions that require the respective heavy modules

import os
import sys
from ngs_pipeline import cerr, cexit, arg_parser


def init_argparser():
    p = arg_parser("convert GFF3 files to BED files")

    p.add_argument(
        "-o", "--outfile", default="-", help="Output filename (default: stdout)"
    )
    p.add_argument(
        "-f",
        "--field",
        default=[],
        action="append",
        help="Field to be included in the BED file. "
        "Multiple fields can be specified. "
        "Default: all fields",
    )
    p.add_argument("infile", help="source GFF3 file")

    return p


def convert_gff3_to_bed(args):

    from ngs_pipeline import gff3utils
    import pandas as pd

    _, df = gff3utils.gff3_to_df(args.infile)

    # filter the fields
    if any(args.field):
        df = df.loc[df["type"].isin(args.field), :]

    # check for overlapping regions, and if exists, merge them
    last_end = None
    last_start = None
    last_seqid = None
    regions = []
    for i, r in df.iterrows():

        if last_end is None:
            last_start = r["start"]
            last_end = r["end"]
            last_seqid = r["seqid"]
            continue

        if r["seqid"] != last_seqid:
            regions.append((last_seqid, last_start, last_end))
            last_start = r["start"]
            last_end = r["end"]
            last_seqid = r["seqid"]
            continue

        # if the current region overlaps with the previous region, merge them
        if r["start"] < last_end:
            last_end = r["end"]
        else:
            regions.append((last_seqid, last_start, last_end))
            last_start = r["start"]
            last_end = r["end"]

    regions.append((last_seqid, last_start, last_end))

    bed_df = pd.DataFrame(regions, columns=["chrom", "start", "end"])
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)

    if args.outfile == "-":
        outfile = sys.stdout
    else:
        outfile = open(args.outfile, "w")
    bed_df.to_csv(outfile, sep="\t", index=False, header=False)


def main(args):
    convert_gff3_to_bed(args)


# EOF
