# initialize.py - ngs-pipeline command
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

# to improve the responsiveness during bash autocomplete, do not import heavy
# modules (such as numpy, pandas, etc) here, but instead import them within the
# functions that require the respective heavy modules

import os
from pathlib import Path
from ngs_pipeline import cexit, cerr, arg_parser


def init_argparser():
    p = arg_parser('setting up reference sequence')

    p.add_argument('-f', '--fasta', required=True,
                   help='reference sequence')
    p.add_argument('-l', '--length', type=int, default=-1,
                   help='maximum genome length to be called by individual '
                   'joint variant caller')
    p.add_argument('-k', '--key', default='regions',
                   choices=['regions', 'contaminant_regions'],
                   help='dictionary key to be used in the YAML file')
    p.add_argument('-n', '--nolength', default=False, action='store_true',
                   help='no need to calculate the length of each sequences')

    p.add_argument('-o', '--outfile', required=True,
                   help='YAML-formated output filename')
    p.add_argument('-p', '--pattern', default=None,
                   help='pattern to match')
    p.add_argument('--outfasta', default='')

    return p


def get_labels_from_fasta(infile, count_length=True):
    """ return list of [(label, length)] """

    sequences = {}
    with open(infile) as f_in:

        curr_label = None
        curr_len = 0
        for line in f_in:
            line = line.strip()
            if line.startswith('>'):
                if curr_label:
                    sequences[curr_label] = curr_len
                curr_label = line[1:].strip().split()[0]
                curr_len = 0
                continue
            if count_length:
                curr_len += len(line.replace(' ',''))

        if curr_label:
            sequences[curr_label] = curr_len

    return sequences


def filter_sequences(infile, outfile, pattern, count_length=False):
    """ return list of [(label, length)] """

    import re

    sequences = {}
    regex = re.compile(pattern)
    with open(infile) as f_in, open(outfile, 'w') as f_out:

        curr_label = None
        curr_len = 0
        save_flag = False
        for line in f_in:
            stripped_line = line.strip()
            if stripped_line.startswith('>'):
                if curr_label:
                    if save_flag:
                        sequences[curr_label] = curr_len
                curr_label = stripped_line[1:].strip().split()[0]
                curr_len = 0
                save_flag = True if regex.match(curr_label) else False
                if save_flag:
                    f_out.write(line)
                continue

            if count_length:
                curr_len += len(stripped_line.replace(' ',''))

            if save_flag:
                f_out.write(line)

        if curr_label:
            if save_flag:
                sequences[curr_label] = curr_len

    return sequences


def generate_ranges(total_length, avg_length):
    rv = round(total_length / avg_length / 1.1)
    if rv == 0:
        rv = 1
    actual_avg_length = round(total_length / rv) + 1
    return list(range(1, total_length, actual_avg_length))[1:] + [total_length]


def partition_regions(sequences, length):
    partitions = {}
    for r, l in sequences.items():
        partitions[r] = generate_ranges(l, length)

    return partitions


def setup_references(args):

    import yaml
    
    partitioned = False
    if args.pattern and args.outfasta:
        seq_labels = filter_sequences(args.fasta, args.outfasta,
                                      args.pattern, True)
    else:
        seq_labels = get_labels_from_fasta(args.fasta, not args.nolength)

    if args.length > 0:
        seq_labels = partition_regions(seq_labels, args.length)
        partitioned = True

    if args.nolength:
        seq_labels = list(seq_labels.keys())

    yaml.dump({args.key: seq_labels}, open(args.outfile, 'w'),
              default_flow_style=None if partitioned else False)


def main(args):
    setup_references(args)

# EOF
