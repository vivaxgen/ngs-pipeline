__copyright__ = '''
generate_variant_report.py - ngs-pipeline command line
[https://github.com/vivaxgen/ngs-pipeline]

(c) 2023 Hidayat Trimarsanto <trimarsanto@gmail.com>
(c) 2024 Ludwig Hoon <ldwgkshoon@gmail.com>

All right reserved.
This software is licensed under MIT license.
Please read the README.txt of this software.
'''

# to improve the responsiveness during bash autocomplete, do not import heavy
# modules (such as numpy, pandas, etc) here, but instead import them within the
# functions that requires respective heavy modules

import sys
import os
from ngs_pipeline import cerr, arg_parser


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
    p.add_argument('--clair3_gvcf', action='store_true',
                   help='input file is in Clair3 gvcf format')
    p.add_argument('--flag_failed_variant', action='store_true',
                   help='Output all variant, including those with depth < mindepth marked with (*) and qual < minqual marked with (^)')
    return p


def generate_variant_report(args):

    # import heavy modules here if required
    from pathlib import Path
    import pandas as pd
    import sys
    from cyvcf2 import VCF
    import IPython

    cerr(f'[Reading variant information from {args.infofile}]')

    info_df = pd.read_table(args.infofile)
    info_df["var_name"] = info_df.Name + " (" + info_df["CHROM"] + ":" + info_df["POS"].astype(str) + ":" + info_df["Change"] + ")"
    info_df.set_index(['CHROM', 'POS'], inplace=True)
    info_df.sort_index(inplace=True)
    

    cerr(f'[Generating variant report from {args.infile}]')

    variants = []
    alleles = []

    vcf = VCF(args.infile, gts012=True)
    sample = vcf.samples[0]
    is_clair3_gvcf = args.clair3_gvcf
    for v in vcf:

        try:

            info_rows = info_df.loc[(v.CHROM, v.POS)] # could be multiple mutations at the same position
            
            var_name = list(info_rows.var_name)
            var_change = list(info_rows.Change)
            variants.extend(var_name)
            # Get the interested variant (which could be ref)

            if is_clair3_gvcf:
                # if GT != 0/0, assert that the ALT is not <NON_REF>
                ALT = max(v.genotypes[0][0:2])
                # Sanity check, if failed, something must have gone wrong
                if ALT > 0:
                    assert([v.REF, *v.ALT][ALT] != "<NON_REF>")

            if 'DP' in v.FORMAT:
                allele_depth = v.format('DP')[0][0]
            elif 'MIN_DP' in v.FORMAT:
                allele_depth = v.format('MIN_DP')[0][0]
            else:
                try:
                    allele_depth = v.INFO['DP']
                except KeyError:
                    allele_depth = 0

            # check if depth is sufficent
            depth_fail = allele_depth < args.mindepth
            
            # if no reads at all, mark as unknown
            # OR if there is no need to flag failed variant AND depth is too low
            if ((not args.flag_failed_variant) and depth_fail) or allele_depth == 0:
                alleles.extend(['?' for _ in var_name])
                continue

            # if low variant quality, then mark as ?, since the interested variant could be the ref as well
            # QUAL = 1 - P(locus is homozygous given the data)
            # only check this if the variant is not ref
            qual_fail = v.QUAL < args.min_var_qual

            gt_changes = gts_to_gt_changes(v.genotypes[0][:-1], v.REF, v.ALT) # can GT be 1/2? e.g., 154535277	T>C and 154535277	T>A at the same time?

            interested_variants = gt_changes & set(var_change)
            
            if gt_changes == {'?'}:
                alleles.extend(['?' for _ in var_name])
                continue

            # No interested variant
            if len(interested_variants) == 0:
                # No interested variant, but depth too low to be confident, will only come to here if --flag_failed_variant
                if depth_fail:
                    alleles.extend(['-*' for _ in var_name])
                    continue
                alleles.extend(['-' for _ in var_name])
                continue
            
            non_ref_variants_index = [i for i, change in enumerate(var_change) if not change.endswith("=")]
            # has interested variant
            if len(interested_variants) > 0:
                # hom
                if len(gt_changes) == 1:
                    temp = ["+" if change in interested_variants else "-" for change in var_change ]

                # het
                elif len(gt_changes) > 1:
                    temp = ["-/+" if change in interested_variants else "-" for change in var_change ]

                # Flag variant that failed mindepth and minqual parameters
                if args.flag_failed_variant:
                    if depth_fail:
                        temp = [var + "*" for var in temp]
                    if qual_fail:
                        temp = [var + "^" if i in non_ref_variants_index else var for i, var in enumerate(temp) ]
                # flag_failed_variant is not set, low depth variant is handled above, these variants will be lowqual
                # Change "+" to "-" if the variant quality is low
                elif qual_fail:
                    temp = ["-" if i in non_ref_variants_index else var for i, var in enumerate(temp) ]

                alleles.extend(temp)
                continue

            # Didn't catch any condition, Sanity check
            alleles.append('!')
            # IPython.embed()
        
        except KeyError:
            cerr(f'[WARNING] is_clair_gvcf:{is_clair3_gvcf} - {v.CHROM}:{v.POS} '
                 f'REF: {v.REF}, ALT: {v.ALT} is not found in the infofile, skipping')


    # Append those variants that are not present in the VCF
    variant_not_in_vcf = [a for a in list(info_df.var_name) if not a in variants]
    variants.extend(variant_not_in_vcf)
    alleles.extend(['?' for _ in variant_not_in_vcf])

    report_df = pd.DataFrame([[sample] + alleles], columns=['SAMPLE'] + variants)
    ordered_columns = ['SAMPLE']
    ordered_columns.extend(info_df.var_name.tolist())
    report_df = report_df[ordered_columns]
    report_df.to_csv(args.outfile, index=False, sep='\t')



def gts_to_gt_changes(gts, ref, alts):
    changes = []
    for gt in gts:
        if gt == -1:
            changes.append("?")
        elif gt == 0:
            changes.append(ref + "=")
        else:
            changes.append(ref + ">" + alts[gt-1])
    return set(changes)

def main(args):
    generate_variant_report(args)

# EOF
