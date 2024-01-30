#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

__copyright__ = '''
run_amtofastq.py - ngs-pipeline command line
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

from ngs_pipeline import cerr, cexit, run_main, arg_parser


# usage: run_amtofastq.py

def init_argparser():
    p = arg_parser(desc='convert alignment map file (SAM/BAM/CRAM) to fastq files')
    p.add_argument('-j', type=int, default=72)
    p.add_argument('--dryrun', default=False, action='store_true')
    p.add_argument('--showcmds', default=False, action='store_true')
    p.add_argument('--unlock', default=False, action='store_true')
    p.add_argument('--rerun', default=False, action='store_true')
    p.add_argument('--target', default='all')
    p.add_argument('--srcdir', default='.',
                   help="Set source directory, where all the source files reside")
    p.add_argument('-s', default=None,
                   help='Sample code, to be used for prefix of fastq filenames')
    p.add_argument('-i', default=None,
                   help='Alignment map file to be splitted to fastq files')
    p.add_argument('infile', default=None, nargs='?',
                   help='A tab-separated value with header SAMPLE and SOURCE')
    return p


def run_amtofastq(args):

    import pathlib
    import snakemake
    from pathlib import Path

    # set SOURCE data

    if args.s:
        if not args.i:
            cexit('Please provide -i if you use -s')
        sources = {args.s: (Path(args.srcdir) / args.i).as_posix()}

    elif args.infile:
        import pandas as pd
        df = pd.read_table(args.infile)
        sources = {}
        for idx, row in df.iterrows():
            sources[row['SAMPLE']] = (Path(args.srcdir) / row['SOURCE']).as_posix()

    else:
        cexit('Please either provide INFILE or use -s & -i')

    # run smk

    snakemake.snakemake(
        pathlib.Path(os.environ['NGS_PIPELINE_BASE']) / 'rules' / 'amtofastq.smk',
        config=dict(SOURCES=sources),
        dryrun=args.dryrun,
        printshellcmds=args.showcmds,
        unlock=args.unlock,
        force_incomplete=args.rerun,
        cores=args.j,
        cluster=os.environ.get('JOBCMD', ''),
        cluster_cancel="scancel",
        targets=[args.target],
    )


if __name__ == '__main__':
    run_main(init_argparser, run_amtofastq)

# EOF
