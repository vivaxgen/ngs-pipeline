__copyright__ = '''
generate_manifest.py - ngs-pipeline command line
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

    m = p.add_mutually_exclusive_group()
    m.add_argument('-s', '--single', default=False, action='store_true',
                   help='fastq files are single (non-paired) such as ONT reads')
    m.add_argument('-p', '--paired', default=False, action='store_true',
                   help='fastq files are paired such as Illumina paired-end reads')

    p.add_argument('-u', '--underscore', type=int, default=0,
                   help='number of consecutive underscore to be stripped from'
                   'filenames to form sample code, counted in reverse')

    p.add_argument('--pause', type=int, default=0,
                   help='pause (in seconds) before back to shell prompt, '
                   'useful for automatic or batch processing so users can '
                   'double check the sample names')
    p.add_argument('--ask-confirmation', default=False, action='store_true',
                   help='ask confirmation to continue saving to output file')

    p.add_argument('-i', '--initial-manifest', default=None, type=str,
                   help='optional initial manifest file to be appended')

    p.add_argument('infiles', nargs='*')
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

    initial_df = None
    if args.initial_manifest:
        cerr(f'[Reading initial manifest file: {args.initial_manifest}]')
        initial_df = pd.read_table(args.initial_manifest)
        if not (('SAMPLE' in initial_df.columns) and ('FASTQ' in initial_df.columns)):
            cexit('[Initial manifest file does not have SAMPLE and/or FASTQ columns]')

    cerr(f'[Receiving {len(args.infiles)} source files]')
    if (initial_df is None) and (len(args.infiles) == 0):
        cexit('[No source files as input, please check your input files]')

    if args.single:
        mode = fileutils.ReadMode.SINGLETON
    elif args.paired:
        mode = fileutils.ReadMode.PAIRED_END
    else:
        mode = None

    read_files = fileutils.ReadFileDict(args.infiles, args.underscore, mode)

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
    if initial_df is not None:
        df = pd.concat([initial_df, df])

    # print at least 5 rows
    cerr('[Showing snippets of sample(s):]')
    cerr(str(df))
    cerr('If the sample code is not correct, try to use option --underscore')

    if args.pause > 0:
        import time
        cerr(f'[Pausing for {args.pause} second(s) for manual inspection]')
        cerr('[Press CTRL-C to abort saving the output file]')
        time.sleep(args.pause)
    
    if args.ask_confirmation:
        resp = input(f'Continue saving to {args.outfile} [y/n]: ')
        if resp.lower().strip()[0] != 'y':
            return

    df.to_csv(args.outfile, sep='\t', index=False)
    cerr(f'[Writing {len(df)} sample manifest to {args.outfile}]')


def main(args):
    generate_manifest(args)

# EOF
