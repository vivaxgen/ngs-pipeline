__copyright__ = '''
gather_stats.py - ngs-pipeline command line
[https://github.com/vivaxgen/ngs-pipeline]

(c) 2022 Hidayat Trimarsanto <trimarsanto@gmail.com>

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
    p = arg_parser('gather stats')
    p.add_argument('-s', '--statfile', default='logs/stats.tsv',
                   help='stat file to be gathered [logs/stats.tsv]')
    p.add_argument('-o', '--outfile', required=True,
                   help='output filename, eg. my-stats.tsv')
    p.add_argument('indirs', nargs='+')
    return p


def gather_stats(args):

    # import heavy modules here if required
    from pathlib import Path
    import pandas as pd

    cerr('Gathering stats from pipeline results')

    dfs = []

    for indir in args.indirs:
        for statfile in Path(indir).glob('*/' + args.statfile):
            dfs.append(pd.read_table(statfile, sep='\t'))

    # concat all
    cerr(f'[Gathering stats file from {len(dfs)} directories]')
    all_df = pd.concat(dfs)

    all_df.to_csv(args.outfile, index=False, sep='\t')


def main(args):
    gather_stats(args)

# EOF
