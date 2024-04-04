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

    p.add_argument('-o', '--outfile', required=True,
                   help='YAML-formated output filename')
    return p


def get_labels_from_fasta(infile):
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
            curr_len += len(line.replace(' ',''))

        if curr_label:
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
    seq_labels = get_labels_from_fasta(args.fasta)
    if args.length > 0:
        seq_labels = partition_regions(seq_labels, args.length)
        partitioned = True
    yaml.dump({'regions': seq_labels}, open(args.outfile, 'w'),
              default_flow_style=None if partitioned else False)


def main(args):
    setup_references(args)

# EOF
