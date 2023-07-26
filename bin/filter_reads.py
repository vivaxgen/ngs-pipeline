#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

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

# check that we have NGS_PIPELINE_BASE environemt
if 'NGS_PIPELINE_BASE' not in os.environ:
    print('ERROR: please set proper shell enviroment by sourcing activate.sh',
          file=sys.stderr)
    sys.exit(1)

from ngsutils import cerr, run_main, arg_parser


def init_argparser():
    p = arg_parser('filter reads that map to certain chromosome from stdin to stdout')

    p.add_argument('regions', nargs='+')
    return p


def filter_reads(args):

    # import heavy modules here if required
    import pysam

    regions = set(args.regions)
    cerr('Filtering for {len(regions)} region(s)')

    inbam = pysam.AlignmentFile('-', 'r')
    outbam = pysam.AlignmentFile('-', 'w', template=inbam)

    for aln in inbam:
        if aln.reference_name not in regions:
            continue
        outbam.write(aln)


if __name__ == '__main__':
    run_main(init_argparser, filter_reads)

# EOF
