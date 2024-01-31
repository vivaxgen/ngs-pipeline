__copyright__ = '''
convert_panelinfo_to_bed.py - ngs-pipeline command line
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
    p = arg_parser('convert panel info file to BED')
    p.add_argument('-o', '--outfile', required=True,
                   help='output BED filename, eg. my-position.bed')
    p.add_argument('infile',
                   help='input panel info file, in TSV format')
    return p


def convert_panelinfo_to_bed(args):

    # import heavy modules here if required
    from pathlib import Path
    import pandas as pd

    cerr(f'[Reading variant information from {args.infile}]')

    info_df = pd.read_table(args.infile)
    bed_df = pd.DataFrame(
        dict(
            CHROM=info_df.CHROM,
            START=info_df.POS - 1,
            END=info_df.POS,
            ID=info_df.dbSNP
        )
    )

    bed_df.sort_values(by=['CHROM', 'START'], inplace=True)
    bed_df.to_csv(args.outfile, index=False, header=False, sep='\t')


def main(args):
    convert_panelinfo_to_bed(args)

# EOF
