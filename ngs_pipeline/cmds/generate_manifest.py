__copyright__ = '''
greet.py - ngs-pipeline command line
[https://github.com/vivaxgen/ngs-pipeline]

(c) 2022-2023 Hidayat Trimarsanto <trimarsanto@gmail.com>

All right reserved.
This software is licensed under MIT license.
Please read the README.txt of this software.
'''

# to improve the responsiveness during bash autocomplete, do not import heavy
# modules (such as numpy, pandas, etc) here, but instead import them within the
# functions that requires respective heavy modules

import sys
import os
import pathlib
from ngs_pipeline import cexit, cerr, arg_parser


def init_argparser():
    p = arg_parser('generate sample manifest file')
    p.add_argument('-o', '--outfile', default='outfile.tsv')
    p.add_argument('-s', '--single', default=False, action='store_true',
                   help='fastq files are single (non-paired) such as ONT reads')
    p.add_argument('-u', '--underline', type=int, default=0,
                   help='no of consecutive underlines to be stripped from filenames '
                   'to form sample code, counted in reverse')
    p.add_argument('infiles', nargs='+')
    return p


def get_sample_name(filename, underline=0):
    filename = pathlib.Path(filename)
    if underline != 0:
        return filename.name.rsplit('_', underline)[0]
    return filename.name.removesuffix('.fastq.gz')


def generate_manifest(args):

    # import heavy modules here if required
    import pandas as pd
    from ngs_pipeline import fileutils

    cerr(f'[Receiving {len(args.infiles)} source files]')

    read_files = fileutils.ReadFileDict(args.infiles, args.underline)

    if any(read_files.err_files):
        cexit('ERROR: invalid input files: \n' +
            '\n'.join(f'  {errmsg}' for errmsg in read_files.err_files))

    # for each sample, process manifest line
    sample_series = []
    fastq_series = []
    counter = 0
    for sample in read_files.samples():
        sample_series.append(sample)
        items = [
            ','.join(item) if type(item) == tuple else item for item in read_files[sample]
        ]
        fastq_series.append(';'.join(items))

    df = pd.DataFrame(dict(SAMPLE=sample_series, FASTQ=fastq_series))
    df.to_csv(args.outfile, sep='\t', index=False)

    cerr(f'[Writing {len(df)} sample manifest to {args.outfile}]')


def main(args):
    generate_manifest(args)

# EOF
