__copyright__ = '''
generate_variant_report.py - ngs-pipeline command line
[https://github.com/vivaxgen/ngs-pipeline]

(c) 2023 Hidayat Trimarsanto <trimarsanto@gmail.com>

All right reserved.
This software is licensed under MIT license.
Please read the README.txt of this software.
'''

# to improve the responsiveness during bash autocomplete, do not import heavy
# modules (such as numpy, pandas, etc) here, but instead import them within the
# functions that requires respective heavy modules

import sys
import os
from ngsutils import cerr, arg_parser


def init_argparser():
    p = arg_parser('generate variant reports from VCF')
    p.add_argument('--infofile', required=True,
                   help='a TSV file containing variant information')
    p.add_argument('-o', '--outfile', required=True,
                   help='output filename, eg. my-report.tsv')
    p.add_argument('--mindepth', type=int, default=10,
                   help='minimum depth for calling variants')
    p.add_argument('--min-var-qual', type=float, default=30,
                   help='minimum quality for variants')
    p.add_argument('infile',
                   help='input file in vcf.gz format')
    return p


__done__ = False
def d():
    global __done__
    __done__ = True


def generate_variant_report(args):

    # import heavy modules here if required
    from pathlib import Path
    import pandas as pd
    import sys
    from cyvcf2 import VCF
    import IPython

    cerr(f'[Reading variant information from {args.infofile}]')

    info_df = pd.read_table(args.infofile)
    info_df.set_index(['CHROM', 'POS'], inplace=True)

    cerr(f'[Generating variant report from {args.infile}]')

    variants = []
    alleles = []

    vcf = VCF(args.infile, gts012=True)
    sample = vcf.samples[0]

    for v in vcf:

        if __done__:
            break

        try:
            info_row = info_df.loc[(v.CHROM, v.POS)]

            #IPython.embed()

            variants.append(info_row.Name)

            # check if depth is sufficent
            if v.INFO.get('DP') < args.mindepth:
                alleles.append('?')
                continue

            # check if allele is still ref, eg no variant
            if v.gt_types[0] == 0:
                alleles.append('-')
                continue

            # if low variant quality, then assume ref / no variant
            if v.QUAL < args.min_var_qual:
                alleles.append('-')
                continue

            # from now on, this is either loq qualhomo alternate of hets
            if v.gt_types[0] == 2:
                alleles.append('+')
                continue

            if v.gt_types[0] == 1:
                alleles.append('-/+')
                continue

            alleles.append('!')
            # IPython.embed()
        
        except KeyError:
            raise

    report_df = pd.DataFrame([[sample] + alleles], columns=['SAMPLE'] + variants)
    report_df.to_csv(args.outfile, index=False, sep='\t')


def main(args):
    generate_variant_report(args)

# EOF
