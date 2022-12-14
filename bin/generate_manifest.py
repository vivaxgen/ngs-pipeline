#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

__copyright__ = '''
greet.py - ngs-pipeline command line
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
    p = arg_parser('generate sample manifest file')
    p.add_argument('-o', '--outfile', default='outfile.tsv')
    p.add_argument('-s', '--single', default=False, action='store_true',
                   help='fastq files are single (non-paired) such as ONT reads')
    p.add_argument('-u', '--underline', type=int, default=0)
    p.add_argument('infiles', nargs='+')
    return p


def get_sample_name(filename, underline=0):
    if underline != 0:
        return filename.rsplit('_', underline)[0]


def generate_manifest(args):

    # import heavy modules here if required
    import pandas as pd

    if not args.single:
        if (len(args.infiles) % 2) != 0:
            cexit(f'Error: odd number of infiles ({len(args.infiles)})')

    infiles = sorted(args.infiles)

    samples = {}

    if not args.single:

        infiles_1, infiles_2 = infiles[::2], infiles[1::2]
        for (infile_1, infile_2) in zip(infiles_1, infiles_2):
            prefix_1 = get_sample_name(infile_1, underline=args.underline)
            prefix_2 = get_sample_name(infile_2, underline=args.underline)
            if prefix_1 != prefix_2:
                cexit(f'ERROR: unmatch pair [{prefix_1}] <> [{prefix_2}]')

            if prefix_1 not in samples:
                samples[prefix_1] = [f'{infile_1},{infile_2}']
            else:
                samples[prefix_1].append(f'{infile_1},{infile_2}')

    # write to TSV file
    sample_series = []
    fastq_series = []
    for sample in sorted(samples.keys()):
        sample_series.append(sample)
        fastq_series.append(';'.join(samples[sample]))

    df = pd.DataFrame(dict(SAMPLE=sample_series, FASTQ=fastq_series))
    df.to_csv(args.outfile, sep='\t', index=False)


if __name__ == '__main__':
    run_main(init_argparser, generate_manifest)

# EOF
