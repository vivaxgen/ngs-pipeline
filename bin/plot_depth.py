#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

__copyright__ = '''
plot_depth.py - ngs-pipeline command line
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

from ngsutils import cexit, cerr, run_main, arg_parser


def init_argparser():
    p = arg_parser(desc='calculate depth from output of samtools depth')
    p.add_argument('--outplot', default='depths.pdf',
                   help='output plot filename [depths.pdf]')
    p.add_argument('-o', '--outfile', default='-',
                   help='the name of the output file [stdout]')
    p.add_argument('--chroms', default='',
                   help='chromosome to be plotted')
    p.add_argument('--sort', default=False, action='store_true',
                   help='sort chromosome/region order alphabetically')
    p.add_argument('-i', '--infile', default='-',
                   help='the name of file from samtools depth output [stdin]')
    p.add_argument('-d', '--mindepth', type=int, default=5,
                   help='minimum depth to be considered [5]')
    p.add_argument('sample',
                   help='sample name or code')
    return p


def plot_depth(args):

    import pandas as pd
    import numpy as np
    from matplotlib import pyplot as plt

    chroms = args.chroms.split(',')
    if not any(chroms):
        cexit('Please provide chromosome names to be plotted')
    if args.sort:
        chroms = sorted(chroms)

    cerr(f'Reading file {args.infile}')
    try:
        df = pd.read_table(args.infile, header=None)
        df.rename(columns={0: 'CHROM', 1: 'POS', 2: 'DEPTH'}, inplace=True)
    except pd.errors.EmptyDataError:
        df = pd.DataFrame({'CHROM': [], 'POS': [], 'DEPTH': []})

    cerr('Plotting depths...')
    if chroms == ['*']:
        chroms = list(sorted(df.CHROM.unique()))

    R = len(chroms)
    fig, axs = plt.subplots(R, 1, figsize=(25, 5 * R))
    color = 'cornflowerblue'
    if R == 1:
        axs = [axs]

    for ax, chrom in zip(axs, chroms):
        current_df = df[df.CHROM == chrom]
        ax.fill_between(current_df.POS, current_df.DEPTH, facecolor=color, color=color)
        q995 = np.quantile(current_df.DEPTH, 0.995) if len(current_df) > 0 else 1
        if current_df.DEPTH.max() > q995 * 5:
            ax.set_ylim([0, q995 * 5])
        ax.set_ylabel(chrom, fontsize=18)
        # ax.set_yscale('log')

    fig.suptitle(args.sample, fontsize=26)
    fig.tight_layout()
    fig.subplots_adjust(top=0.97)
    fig.savefig(args.outplot)


if __name__ == '__main__':
    run_main(init_argparser, plot_depth)


# EOF
