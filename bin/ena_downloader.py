#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

__copyright__ = '''
ena_downloader.py - ngs-pipeline command line
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

from ngsutils import cerr, cexit, run_main, arg_parser


# to get URL, fetch JSON from the following:
# https://www.ebi.ac.uk/ena/portal/api/filereport?
# result=read_run&fields=fastq_ftp&format=JSON&accession=ERR484725

def get_urls(ena_accid):

    import requests
    import pathlib

    payload = dict(result='read_run',
                   fields='fastq_ftp',
                   format='JSON',
                   accession=ena_accid,
                   )
    r = requests.get('https://www.ebi.ac.uk/ena/portal/api/filereport',
                     params=payload)

    if r.status_code != 200:
        raise ValueError(f'ENA accession not found: {ena_accid}')

    result_list = r.json()
    urls = result_list[0]['fastq_ftp'].split(';')
    filenames = [pathlib.Path(url).name for url in urls]
    return urls, filenames


def fetch_urls(urls):
    """ fetch the urls in parallel """

    import subprocess

    opts = []
    for url in urls:
        opts.append('-O')
        opts.append(url)

    cmds = ['curl', '-C', '-', '-Z'] + opts
    #print(' '.join(cmds))
    subprocess.call(cmds)


def init_argparser():
    p = arg_parser(desc='download FASTQ files based on ENA accession')
    p.add_argument('--outdir', default='.',
                   help='directory where manifest.tsv and FASTQ files would be written to')
    p.add_argument('--outfile', default='out-manifest.tsv')
    p.add_argument('infile',
                   help='tab-separated file to get sample and ENA accession')
    return p


def ena_downloader(args):

    import pandas as pd
    import pathlib

    # sanity check
    outfile = pathlib.Path(args.outfile)
    if outfile.exists():
        cexit(f'ERROR: output file {args.outfile} is already existed.')

    if ':' in args.infile:
        infile, column_specs = args.infile.split(':')
        sample_column, ena_column = column_specs.split(',')
    else:
        infile = args.infile

    cerr(f'[Parsing input file: {infile}]')
    df = pd.read_table(infile, sep=None, engine='python')

    sample_manifest = []
    fastq_manifest = []

    for idx, r in df.iterrows():

        sample = r[sample_column]
        ena_accs = r[ena_column]
        if not ena_accs:
            continue

        cerr(f'[Fetching for sample: {sample}]')
        url_list = []
        filename_list = []
        for ena_acc in ena_accs.split(','):
            cerr(f'[Fetching urls for accession: {ena_acc}]')
            urls, filenames = get_urls(ena_acc)
            url_list += urls
            filename_list.append(','.join(filenames))

        cerr(f'[Parallel downloading for {len(url_list)} urls]')
        fetch_urls(url_list)

        sample_manifest.append(sample)
        fastq_manifest.append(';'.join(filename_list))

    out_df = pd.DataFrame({'SAMPLE': sample_manifest, 'FASTQ': fastq_manifest})
    out_df.to_csv(args.outfile, sep='\t', index=False)
    cerr(f'[Output written to {args.outfile}]')


if __name__ == '__main__':
    run_main(init_argparser, ena_downloader)

# EOF
