#!/usr/bin/env ngs-pl

import sys
import pathlib
from ngs_pipeline import arg_parser, cerr, cexit


def init_argparser():
    p = arg_parser()
    p.add_argument("-o", "--outfile", required=True)
    p.add_argument('infiles', nargs='+')
    return p


def consolidate_haplotypes(unique_haps):

    positions = set()
    
    # convert haplotypes to list of tuples: [(pos, allele), ...]
    # and identify variant positions
    tupled_haplotypes = {}
    for h, k in unique_haps.items():
        if h == '':
            tupled_haplotypes[k] = []
            continue
        h = h.split('|')
        tupled_hap = [(int(t[0]), t[1]) for t in [a.split(':') for a in h]]
        for pos, _ in tupled_hap:
            positions.add(pos)
        tupled_haplotypes[k] = tupled_hap

    # sort variant positions
    positions = sorted(list(positions))

    import IPython; IPython.embed()


def merge_haplotypes(unique_haps):

    import pandas as pd

    # make a list of dictionary
    dicted_haplotypes = []
    for h, k in unique_haps.items():
        if h == '':
            dicted_haplotypes.append(dict(hap_id=k))
            continue
        h = h.split("|")
        dicted_hap = {int(t[0]): t[1] for t in [a.split(':') for a in h]}
        dicted_hap['hap_id'] = k
        dicted_haplotypes.append(dicted_hap)

    df = pd.DataFrame(dicted_haplotypes)

    # reorder columns
    pos_cols = sorted(filter(lambda x: isinstance(x, int), df.columns))
    str_cols = sorted(filter(lambda x: not isinstance(x, int), df.columns))
    df = df[str_cols + pos_cols]
    df.fillna('.', inplace=True)
    
    # print('|'.join(str(x) for x in df.columns))
    # for idx, r in df.iterrows():
    #     c = [f"{i[1]:6}" for i in r.items()]
    #     print('|'.join(c))

    return df
 

def gather_haplotypes(args):

    import json
    import pandas as pd

    hap_list = []
    unique_haps = {}
    last_uniqid = 0

    dfs = []

    for infile in args.infiles:

        infile = pathlib.Path(infile)
        sample = infile.parent.name

        haplotypes = []

        with open(infile) as f:
            cerr(f"Reading {sample=}")
            for line in f:
                tokens = line.split()
                match len(tokens):
                    case 4:
                        haplotype, ratio, depth, sample = tokens
                    case 3:
                        haplotype = ""
                        ratio, depth, sample = tokens
                    case _:
                        raise ValueError(f"infile line is not formatted correctly")

                depth = int(depth)
                ratio = float(ratio)

                try:
                    hapid = unique_haps[haplotype]
                except:
                    hapid = last_uniqid = last_uniqid + 1
                    unique_haps[haplotype] = hapid

                hap_list.append((sample, hapid, depth, ratio, haplotype))

    merge_haplotypes(unique_haps)

    # import pprint
    # pprint.pprint(unique_haps)
    #pprint.pprint(hap_list)
    df = pd.DataFrame(hap_list, columns=['SAMPLE', 'HAPID', 'DEPTH', 'RATIO', 'HAPLOTYPE'])
    df = df.merge(df.SAMPLE.value_counts().to_frame(), on="SAMPLE")
    #import IPython; IPython.embed()
    df.to_csv(args.outfile, sep='\t', index=False)


def main(args):
    gather_haplotypes(args)

# EOF

