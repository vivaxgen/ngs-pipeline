__copyright__ = '''
greet.py - ngs-pipeline command line
[https://github.com/vivaxgen/ngs-pipeline]

(c) 2023 Hidayat Trimarsanto <trimarsanto@gmail.com>

All right reserved.
This software is licensed under MIT license.
Please read the README.txt of this software.
'''

# to improve the responsiveness during bash autocomplete, do not import heavy
# modules (such as numpy, pandas, etc) here, but instead import them within the
# functions that requires respective heavy modules

import sys
import os

from ngs_pipeline import cerr, arg_parser


def init_argparser():
    p = arg_parser('count FR, RF, FF, RR and trans orientaion of a bam file')
    p.add_argument('--remove_FF', action='store_true', default=False)
    p.add_argument('--remove_RR', action='store_true', default=False)
    p.add_argument('--remove_RF', action='store_true', default=False)
    p.add_argument('--remove_FR', action='store_true', default=False)
    p.add_argument('--remove_trans', action='store_true', default=False)
    p.add_argument('--remove_unmapped', action='store_true', default=False)
    p.add_argument('--remove_secondary', action='store_true', default=False)
    p.add_argument('--remove_supplementary', action='store_true', default=False)
    p.add_argument('-o', '--outfile', default='-',
                   help='output bam file. default is stdout [-]')
    p.add_argument('--outstat', default='stat-orientation.json',
                   help='JSON-formatted statistic output file [stat-orientation]')
    p.add_argument('infile',
                   help='input bam file')
    return p


def filter_reads_orientation(args):

    import pysam
    import json
    from ngs_pipeline import add_pgline, get_mode

    # prepare flags

    remove_unmapped = 1 << 0
    remove_trans = 1 << 1
    remove_RF = 1 << 2
    remove_FF = 1 << 3
    remove_RR = 1 << 4
    remove_secondary = 1 << 5
    remove_supplementary = 1 << 6
    remove_FR = 1 << 7

    flags = 0
    if args.remove_unmapped:
        flags |= remove_unmapped
    if args.remove_trans:
        flags |= remove_trans
    if args.remove_RR:
        flags |= remove_RR
    if args.remove_FF:
        flags |= remove_FF
    if args.remove_RF:
        flags |= remove_RF
    if args.remove_secondary:
        flags |= remove_secondary
    if args.remove_supplementary:
        flags |= remove_supplementary
    if args.remove_FR:
        flags |= remove_FR

    in_aln = pysam.AlignmentFile(args.infile, get_mode(args.infile, 'r'))
    header = add_pgline(in_aln,
                        dict(ID='filter_reads_orientation.py',
                             PN='filter_reads_orientation.py',
                             CL=' '.join(sys.argv),
                             DS='filtering reads based on orientation')
                        )

    out_aln = pysam.AlignmentFile(args.outfile, get_mode(args.outfile, 'w'), header=header)

    count_FR = count_RF = count_RR = count_FF = count_trans = count_unmapped = 0
    count_secondary = count_supplementary = count_ERR = 0

    aln_pair_container = {}

    for aln in in_aln:

        # check for secondary & supplementary first
        if aln.is_secondary:
            count_secondary += 1
            if flags & remove_secondary:
                continue
        elif aln.is_supplementary:
            count_supplementary += 1
            if flags & remove_supplementary:
                continue

        if aln.query_name in aln_pair_container:
            aln1 = aln_pair_container[aln.query_name]
            aln2 = aln
            del aln_pair_container[aln.query_name]
        else:
            aln_pair_container[aln.query_name] = aln
            continue

        if aln1.is_unmapped or aln2.is_unmapped:
            count_unmapped += 1
            if flags & remove_unmapped:
                continue

        elif aln1.reference_name != aln2.reference_name:
            count_trans += 1
            if flags & remove_trans:
                continue

        elif aln1.is_forward and aln2.is_forward:
            # this is FF
            count_FF += 1
            if flags & remove_FF:
                continue

        elif aln1.is_reverse and aln2.is_reverse:
            count_RR += 1
            if flags & remove_RR:
                continue

        elif aln1.is_forward and aln2.is_reverse:
            if aln1.reference_start < aln2.reference_end:
                count_FR += 1
                if flags & remove_FR:
                    continue
            else:
                count_RF += 1
                if flags & remove_RF:
                    continue
        elif aln1.is_reverse and aln2.is_forward:
            if aln1.reference_end > aln2.reference_start:
                count_FR += 1
                if flags & remove_FR:
                    continue
            else:
                count_RF += 1
                if flags & remove_RF:
                    continue
        else:
            import IPython; IPython.embed()
            count_ERR += 1
            continue

        out_aln.write(aln1)
        out_aln.write(aln2)

    non_paired_segments = list(aln_pair_container.keys())
    d = dict(
        FR_orientation=count_FR,
        RF_orientation=count_RF,
        FF_orientation=count_FF,
        RR_orientation=count_RR,
        trans_orientation=count_trans,
        ERR_orientation=count_ERR,
        secondary_alignment=count_secondary,
        supplementary_alignment=count_supplementary,
        non_paired_N=len(non_paired_segments),
        non_paired_segments=non_paired_segments,
    )
    json.dump(d, open(args.outstat, 'w'))


def main(args):
    filter_reads_orientation(args)

# EOF
