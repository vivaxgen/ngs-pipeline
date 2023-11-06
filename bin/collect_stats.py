#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

__copyright__ = '''
.py - ngs-pipeline command line
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
    p = arg_parser('collect stats from various steps')
    p.add_argument('-o', '--outfile', default='-',
                   help='output filename [stdout]')
    p.add_argument('--trimmed', action='append')
    p.add_argument('--mapped', action='append')
    p.add_argument('--final', action='append')
    p.add_argument('--dedup', action='append')
    p.add_argument('--depth', action='append')
    p.add_argument('--mindepth', default=5, type=int,
                   help='min depth to be stats [5]')
    p.add_argument('sample',
                   help='Sample code')
    return p


def parse_read_trimming(infiles):
    """ return (no_of_original_reads, no_of_filtered_reads) """

    # this will read a json file with the following specs:
    # {'original_reads': int, 'filtered_reads': int}

    import json
    original_reads = 0
    filtered_reads = 0

    for infile in infiles:
        with open(infile, 'r') as f_in:
            d = json.load(f_in)
            original_reads += d['original_reads']
            filtered_reads += d['filtered_reads']

    return (original_reads, filtered_reads)


def parse_number(line, func=int):
    # read number after colon
    return func(line.split(':')[-1].strip().split()[0].replace(',',''))


def parse_mapping_stats(infiles):
    # read from output of samtools stats

    mapped = 0
    proper = 0
    insert_size = 0.0
    insert_size_sd = 0.0
    avg_qual = 0.0
    inward = 0
    outward = 0
    other = 0
    trans = 0

    N = len(infiles)

    for infile in infiles:
        cerr(f'[Reading stat file: {infile}]')
        with open(infile, 'r') as f_in:
            for line in f_in:
                if not line.startswith('SN\t'):
                    continue
                if 'reads mapped:' in line:
                    mapped += parse_number(line)
                    continue
                if 'reads properly paired:' in line:
                    proper += parse_number(line)
                    continue
                if 'insert size average:' in line:
                    insert_size += parse_number(line, float)
                    continue
                if 'insert size standard deviation:' in line:
                    insert_size_sd += parse_number(line, float)
                    continue
                if 'average quality:' in line:
                    avg_qual += parse_number(line, float)
                    continue
                if 'inward oriented pairs:' in line:
                    inward += parse_number(line)
                    continue
                if 'outward oriented pairs:' in line:
                    outward += parse_number(line)
                    continue
                if 'pairs with other orientation:' in line:
                    other += parse_number(line)
                    continue
                if 'pairs on different chromosomes:' in line:
                    trans += parse_number(line)
                    continue

    return (
        mapped, proper, insert_size/N, insert_size_sd/N, avg_qual/N,
        inward * 2, outward * 2, other * 2, trans * 2
    )


def calculate_depth(infile, mindepth=5):

    import gzip
    import statistics

    depths = []
    with gzip.open(infile, 'r') as fin:
        cerr(f'[Reading depth file: {infile}]')
        for line in fin:
            d = int(line.strip().split()[-1])
            if d > mindepth:
                depths.append(d)

    total = sum(depths)
    L = len(depths)
    if total == 0:
        average = q1 = 0
    else:
        average = total / len(depths)
        q1 = statistics.quantiles(depths, n=4)[0]

    return (L, average, q1)


def calculate_depths(infiles, mindepth=5):

    import gzip
    import statistics

    depths = {}
    for infile in infiles:
        with gzip.open(infile, 'r') as fin:
            cerr(f'[Reading depth file: {infile}]')
            for line in fin:
                tokens = line.strip().split()
                try:
                    chr_d = depths[tokens[0]]
                except KeyError:
                    chr_d = depths[tokens[0]] = {}

                d = int(tokens[-1])
                try:
                    chr_d[tokens[1]] += d
                except KeyError:
                    chr_d[tokens[1]] = d

    depth_values = []

    for chr_d in depths.values():
        for d in chr_d.values():
            if d >= mindepth:
                depth_values.append(d)

    total = sum(depth_values)
    L = len(depth_values)
    if total == 0:
        average = q1 = 0
    else:
        average = total / L
        q1 = statistics.quantiles(depth_values, n=4)[0]

    return (L, average, q1)


def collect_stats(args):

    # import heavy modules here if required

    initial_reads, trimmed_reads = parse_read_trimming(args.trimmed)
    trimmed_r = trimmed_reads / initial_reads

    mapped_reads, proper_mapped_reads = parse_mapping_stats(args.mapped)[:2]
    mapped_r = (mapped_reads / trimmed_reads) if trimmed_reads > 0 else 0
    proper_r = (proper_mapped_reads / mapped_reads) if mapped_reads > 0 else 0



    _, proper_dedup_reads, insert_size, insert_size_sd, avg_qual, inward, outward, other, trans = parse_mapping_stats(args.final)
    dedup_r = proper_dedup_reads / proper_mapped_reads if proper_mapped_reads > 0 else 0
    final_r = proper_dedup_reads / trimmed_reads if trimmed_reads > 0 else 0
    if proper_dedup_reads > 0:
        inward_r = inward / proper_dedup_reads
        outward_r = outward / proper_dedup_reads
        other_r = other / proper_dedup_reads
        trans_r = trans / proper_dedup_reads
    else:
        inward_r = outward_r = other_r = trans_r = 0

    basepairs, avg, q1 = calculate_depths(args.depth, args.mindepth)

    outfile = sys.stdout if args.outfile == '-' else open(args.outfile, 'w')
    headers = [
        'SAMPLE', 'RAW',
        'TRIMMED', 'TRIMMED_R',
        'MAPPED', 'MAPPED_R',
        'PROPER', 'PROPER_R',
        'DEDUP-FINAL', 'DEDUP-FINAL_R',
        'FINAL_R',
        'INWARD', 'INWARD_R',
        'OUTWARD', 'OUTWARD_R',
        'OTHER', 'OTHER_R',
        'TRANS', 'TRANS_R',
        'INSERT_SIZE', 'INSERT_SIZE_SD', 'AVG_QUAL',
        'BASEPAIRS', 'AVG_DEPTH', 'Q1_DEPTH'
    ]
    outfile.write('%s\n' % '\t'.join(headers))

    outfile.write(
        f'{args.sample}\t'
        f'{initial_reads}\t'
        f'{trimmed_reads}\t{trimmed_r:5.3f}\t'
        f'{mapped_reads}\t{mapped_r:5.3f}\t'
        f'{proper_mapped_reads}\t{proper_r:5.3f}\t'
        f'{proper_dedup_reads}\t{dedup_r:5.3f}\t'
        f'{final_r:5.3f}\t'
        f'{inward}\t{inward_r:5.3f}\t'
        f'{outward}\t{outward_r:5.3f}\t'
        f'{other}\t{other_r:5.3f}\t'
        f'{trans}\t{trans_r:5.3f}\t'
        f'{insert_size}\t{insert_size_sd}\t{avg_qual}\t'
        f'{basepairs}\t{avg:5.2f}\t{q1:5.2f}'
        '\n'
    )


if __name__ == '__main__':
    run_main(init_argparser, collect_stats)


# EOF

