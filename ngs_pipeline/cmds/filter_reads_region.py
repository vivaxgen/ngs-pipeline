__copyright__ = """
greet.py - ngs-pipeline command line
[https://github.com/vivaxgen/ngs-pipeline]

(c) 2023 Hidayat Trimarsanto <trimarsanto@gmail.com>

All right reserved.
This software is licensed under MIT license.
Please read the README.txt of this software.
"""

# to improve the responsiveness during bash autocomplete, do not import heavy
# modules (such as numpy, pandas, etc) here, but instead import them within the
# functions that requires respective heavy modules

import sys
import os

from ngs_pipeline import cerr, arg_parser


def init_argparser():
    p = arg_parser(
        "filter (in/out) reads that map to certain chromosome from stdin to stdout, "
        "including the mate"
    )
    p.add_argument(
        "--remove",
        default=False,
        action="store_true",
        help="removing (instead of keeping) read that mapped to the regions or whose mate "
        "mapped to the regions, keeping others (including unmapped pairs)",
    )
    p.add_argument("--outstat", default="", help="JSON-based stat output file path")
    p.add_argument(
        "-o", "--outfile", default="-", help="output filename, default is stdout [-]"
    )
    p.add_argument("regions", nargs="+")
    return p


def filter_reads(args):

    # import heavy modules here if required
    import pysam
    import json

    regions = set(args.regions)
    cerr("Filtering for {len(regions)} region(s)")

    inbam = pysam.AlignmentFile("-", "r")

    if args.outfile == "-":
        outbam = pysam.AlignmentFile("-", "w", template=inbam)
    else:
        outbam = pysam.AlignmentFile(args.outfile, "wb", template=inbam)

    # we process alignment segments in pairs
    aln_pairs = {}

    in_counter = 0
    out_counter = 0

    if args.remove:

        cerr("Mode: Filter-out")

        for aln in inbam:

            if aln.query_name not in aln_pairs:
                aln_pairs[aln.query_name] = aln
                continue
            aln_pair = aln_pairs[aln.query_name]
            del aln_pairs[aln.query_name]

            # check that both are mapped, otherwise drop
            if not (aln.reference_name and aln_pair.reference_name):
                out_counter += 1
                continue

            # we remove alignment segment pair if any of the pair maps to regions
            if aln.reference_name in regions or aln_pair.reference_name in regions:
                out_counter += 1
                continue

            outbam.write(aln_pair)
            outbam.write(aln)
            in_counter += 1

    else:

        cerr("Mode: Filter-in")

        for aln in inbam:

            cerr(f"aln.query_name: {aln.query_name}")
            if aln.query_name not in aln_pairs:
                aln_pairs[aln.query_name] = aln
                cerr(f"aln_pairs size: {len(aln_pairs)}")
                continue
            aln_pair = aln_pairs[aln.query_name]
            del aln_pairs[aln.query_name]

            # check that both are mapped, otherwise drop
            if not (aln.reference_name and aln_pair.reference_name):
                out_counter += 1
                continue

            # we write alignment segment pair if all of the pair map to regions

            if aln.reference_name in regions and aln_pair.reference_name in regions:
                outbam.write(aln_pair)
                outbam.write(aln)
                in_counter += 1
                continue

            out_counter += 1

    if args.outstat:
        d = dict(
            mode="remove" if args.remove else "write",
            filter_in=in_counter,
            filter_out=out_counter,
            error_pairs=len(aln_pairs),
        )
        json.dump(d, open(args.outstat, "w"))


def filter_reads_region(args):

    # import heavy modules here if required
    import pysam
    import json
    from ngs_pipeline import get_mode, add_pgline

    regions = set(args.regions)
    cerr(f"Filtering for {len(regions)} region(s)")

    # note that the input file has to be pre-processed to fix mate information
    # either using samtools fixmate, picard or any other
    inbam = pysam.AlignmentFile("-", "r")

    header = add_pgline(
        inbam,
        dict(
            ID="filter_reads_region.py",
            PN="filter_reads_region.py",
            CL=" ".join(sys.argv),
            DS="filtering reads based on regions",
        ),
    )
    outbam = pysam.AlignmentFile(
        args.outfile, get_mode(args.outfile, "w"), header=header
    )

    in_counter = 0
    out_counter = 0

    if args.remove:

        for aln in inbam:

            if aln.reference_name in regions or aln.next_reference_name in regions:
                out_counter += 1
                continue

            outbam.write(aln)
            in_counter += 1

    else:

        for aln in inbam:

            # we write alignment segment pair if all of the pair map to regions

            if aln.reference_name in regions and aln.next_reference_name in regions:
                outbam.write(aln)
                in_counter += 1
                continue

            out_counter += 1

    if args.outstat:
        d = dict(
            mode="remove" if args.remove else "write",
            filter_in=in_counter,
            filter_out=out_counter,
        )
        json.dump(d, open(args.outstat, "w"))


def main(args):
    filter_reads_region(args)


# EOF
