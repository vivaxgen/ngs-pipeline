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
    p.add_argument('-o', '--outfile', default=None,
                   help='the name of the output file (eg. keywords.txt)')
    p.add_argument('-c', '--configfile', default=None,
                   help='configuration file to match against')
    p.add_argument('indirs', nargs='+')

    return p


def scan_config_keywords(args):

    import yaml

    keywords = []
    for indir in args.indirs:
        keywords += snakeutils.scan_for_config_keywords(indir)

    if args.outfile:
        with open(args.outfile, 'w') as f_out:
            keywords = sorted(keywords)
            f_out.write('\n'.join(keywords))

    if args.configfile:
        keywords = set(keywords)
        with open(args.configfile) as f_in:
            configs = yaml.safe_load(f_in)
            in_configs = set(configs.keys())
            unused = keywords - in_configs
            unknown = in_configs - keywords

        if any(unused):
            cerr('Unused configuration keys:')
            cerr('\n'.join(sorted(unused)))
            cerr('')

        if any(unknown):
            cerr('Unknown configuration keys:')
            cerr('\n'.join(sorted(unknown)))
            cerr('')


def main(args):
    scan_config_keywords(args)

# EOF
