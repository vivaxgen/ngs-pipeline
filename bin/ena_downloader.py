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
#
# in case fastq_ftp returns empty string, use this instead:
# result=read_run&fields=submitted_ftp&format=JSON&accession=ERR2299660
#
# submitted_ftp is usually a cram or bam file, which then needs to be splitted:
# samtools collate -u -f -O file_cram.cram |
# samtools fastq -1 read_1.fastq.gz -2 read_2.fastq.gz -0 /dev/null -s /dev/null -n
#
# for API docs, see https://www.ebi.ac.uk/ena/portal/api/doc?format=pdf
#

def get_urls(ena_accid, query):

    import requests
    import pathlib
    from urllib.parse import quote

    payload = dict(result='read_run',
                   fields=query,
                   format='JSON',
                   accession=ena_accid,
                   )
    r = requests.get('https://www.ebi.ac.uk/ena/portal/api/filereport',
                     params=payload)

    if r.status_code != 200:
        raise ValueError(f'ENA accession not found: {ena_accid}.')

    result_resp = r.json()
    query_resp = result_resp[0][query]
    if not query_resp:
        raise ValueError(f'ENA accession does not have {query} entry.')

    urls = query_resp.split(';')

    filenames = [pathlib.Path(url).name for url in urls]
    print(urls, filenames)
    return [quote(url) for url in urls], filenames


def fetch_urls(urls, filenames, outdir, socks5=None):
    """ fetch the urls in parallel """

    import subprocess
    import urllib.parse

    opts = []
    for url, filename in zip(urls, filenames):
        opts.append('-o')
        opts.append(filename)
        opts.append(url)

    cmds = ['curl']
    if socks5:
        cmds += ['--socks5-hostname', socks5]

    cmds += ['-C', '-', '-Z']
    cmds += opts

    print(' '.join(cmds))
    return subprocess.call(cmds, cwd=outdir)


def init_argparser():
    p = arg_parser(desc='download FASTQ files based on ENA accession')
    p.add_argument('--socks5', default='',
                   help='socks5 hostname:port to be used as ftp proxy')
    p.add_argument('--outdir', default='.',
                   help='directory where manifest.tsv and FASTQ files would be written to, '
                   'default to current directory')
    p.add_argument('--outfile', default='out-manifest.tsv',
                   help='default to out-manifest.tsv')
    p.add_argument('--outurls', default='url-manifest.txt',
                   help='ftp list file that can be used by eg. wget, default to url-manifest.txt')
    p.add_argument('--outerr', default='err-log.txt',
                   help='default to err-log.txt')
    p.add_argument('--wait', default=0, type=int,
                   help='waiting time (seconds) between requests to lighten the remote server')
    p.add_argument('infile',
                   help='tab-separated file to get sample and ENA accession, use colon to '
                   'notation to indicate the columns to be used, eg: FILE.TSV:SAMPLE,ENA')
    return p


def ena_downloader(args):

    import pandas as pd
    import pathlib
    import time

    # put all output files to out directory
    args.outfile = f'{args.outdir}/{args.outfile}'
    args.outurls = f'{args.outdir}/{args.outurls}'
    args.outerr = f'{args.outdir}/{args.outerr}'

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

    outurls  = open(args.outurls, 'w')
    outerr = open(args.outerr, 'w')

    N = len(df)
    cerr(f'[Start downloading for {N} sample(s)]')

    for idx, r in df.iterrows():

        sample = r[sample_column]
        ena_accs = r[ena_column]
        if not ena_accs or type(ena_accs) != str:
            continue

        cerr(f'[{idx+1} of {N}: Fetching for sample: {sample}]')
        url_list = []
        filename_list = []
        for ena_acc in ena_accs.split(','):
            cerr(f'[Fetching urls for accession: {ena_acc}]')
            ena_acc = ena_acc.strip()
            ok = False
            # we try fastq files first, if no fastq then  whatever submitted to ENA
            for q in ['fastq_ftp', 'submitted_ftp']:
                try:
                    urls, filenames = get_urls(ena_acc, q)
                    ok = True
                    break
                except ValueError as exc:
                    outerr.write(f'WARN: {exc}\n')
                    continue
            url_list += urls
            outurls.write('\n'.join(urls))
            outurls.write('\n')
            filename_list.append(filenames)

        cerr(f'[Parallel downloading for {len(url_list)} urls]')
        ret = fetch_urls(url_list, sum(filename_list, []), args.outdir, args.socks5)
        if ret != 0:
            cerr(f'[WARN: non-zero return code from curl]')
            outerr.write(f'WARN: non-zero return code from curl for URLS: {",".join(url_list)}')

        sample_manifest.append(sample)
        fastq_manifest.append(';'.join(','.join(paired_file) for paired_file in filename_list))

        if args.wait > 0:
            time.sleep(args.wait)

    out_df = pd.DataFrame({'SAMPLE': sample_manifest, 'FASTQ': fastq_manifest})
    out_df.to_csv(args.outfile, sep='\t', index=False)
    cerr(f'[Manifest written to {args.outfile}]')
    outurls.close()
    outerr.close()
    cerr(f'[URL list written to {args.outurls}]')


if __name__ == '__main__':
    run_main(init_argparser, ena_downloader)

# EOF
