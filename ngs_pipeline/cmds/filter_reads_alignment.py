# filter_reads_alignment.py - ngs-pipeline command line
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(c) 2025, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

# to improve the responsiveness during bash autocomplete, do not import heavy
# modules (such as numpy, pandas, etc) here, but instead import them within the
# functions that require the respective heavy modules

import sys
from ngs_pipeline import cout, cerr, arg_parser, check_NGS_PIPELINE_BASE

# this tool filters reads from BAM file based on alignment-specific criteria
# suitable mostly for long reads (ONT, PacBio, etc)


def init_argparser():

    p = arg_parser("filter reads from BAM based on alignment criteria")
    p.add_argument(
        "--min-matching-length",
        type=int,
        default=-1,
        help="filter reads based on minimum matching length (ie total number of M in CIGAR), default is -1 (no limit)",
    )
    p.add_argument(
        "-o", "--outfile", default="-", help="output bam file. default is stdout [-]"
    )
    p.add_argument("infile", help="input bam file")
    return p


def filter_reads_alignment(args):

    import pysam
    from ngs_pipeline import add_pgline, get_mode

    in_aln = pysam.AlignmentFile(args.infile, get_mode(args.infile, "r"))
    header = add_pgline(
        in_aln,
        dict(
            ID="filter_reads_alignment.py",
            PN="filter_reads_alignment.py",
            CL=" ".join(sys.argv),
            DS="filtering reads based on alignment criteria",
        ),
    )

    out_aln = pysam.AlignmentFile(
        args.outfile, get_mode(args.outfile, "w"), header=header
    )

    for aln in in_aln:
        if args.min_matching_length > 0:
            # get total matching length from CIGAR
            if aln.get_cigar_stats()[0][0] < args.min_matching_length:
                continue

        out_aln.write(aln)

    cerr(f"output written to {args.outfile}")


def main(args):
    filter_reads_alignment(args)


# EOF
