# scan_config_keywords.py - ngs-pipeline command
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

# to improve the responsiveness during bash autocomplete, do not import heavy
# modules (such as numpy, pandas, etc) here, but instead import them within the
# functions that require the respective heavy modules


import os
from ngs_pipeline import (cerr, cexit, arg_parser,
                          check_NGSENV_BASEDIR, check_NGS_PIPELINE_BASE,
                          get_snakefile_path, setup_config, snakeutils)

def init_argparser():
    p = arg_parser(desc='scan configuration keywords on snakefiles')
    p.add_argument('-o', '--outfile', default='keywords.txt',
                   help='the name of the output file')
    p.add_argument('indirs', nargs='+')

    return p


def scan_config_keywords(args):

    keywords = []
    for indir in args.indirs:
        keywords += snakeutils.scan_for_config_keywords(indir)

    print(keywords)


def main(args):
    scan_config_keywords(args)

# EOF
