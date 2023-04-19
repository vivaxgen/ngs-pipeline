#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

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

# check that we have NGS_PIPELINE_BASE environemt
if 'NGS_PIPELINE_BASE' not in os.environ:
    print('ERROR: please set proper shell enviroment by sourcing activate.sh',
          file=sys.stderr)
    sys.exit(1)

from ngsutils import cerr, run_main, arg_parser


def init_argparser():
    p = arg_parser('gather stats')
    p.add_argument('-s', '--statfile', default='logs/stats.tsv',
                   help='stat file to be gathered [logs/stats.tsv]')
    p.add_argument('-o', '--outfile',
                   help='output filename')
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


if __name__ == '__main__':
    run_main(init_argparser, gather_stats)

# EOF
