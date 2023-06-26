#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

__copyright__ = '''
count_depth.py - ngs-pipeline command line
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

import sys
import gzip


def init_argparser():
    p = arg_parser(desc='calculate depth from output of samtools depth')
    p.add_argument('-o', '--outfile', default='-',
                   help='the name of the output file [stdout]')
    p.add_argument('-i', '--infile', default='-',
                   help='the name of file from samtools depth output [stdin]')
    p.add_argument('-d', '--mindepth', type=int, default=5,
                   help='minimum depth to be considered [5]')
    p.add_argument('sample',
                   help='sample name or code')
    return p


def count_depth(args):

    depths = []
    if args.infile == '-':
        f_in = sys.stdin
    else:
        f_in = gzip.open(args.infile, 'r')

    for line in f_in:
        d = int(line.strip().split()[-1])
        if d >= args.mindepth:
            depths.append(d)

    total = sum(depths)
    if total == 0:
        avg = 0
    else:
        avg = total / len(depths)

    if args.outfile == '-':
        f_out = sys.stdout
    else:
        f_out = open(args.outfile, 'w')
    print(f'{args.sample} {avg} {len(depths)}', file=f_out)


if __name__ == '__main__':
    run_main(init_argparser, count_depth)


# EOF

