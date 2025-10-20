# showcmds.py - ngs-pipeline command line
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(c) 2023, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

# to improve the responsiveness during bash autocomplete, do not import heavy
# modules (such as numpy, pandas, etc) here, but instead import them within the
# functions that require the respective heavy modules

from ngs_pipeline import cout, cerr, arg_parser
from ngs_pipeline import subcommands


def init_argparser():

    p = arg_parser("show all available commands")
    return p


def main(args):

    subcommands.get_subcommands().show_commands()


# EOF
