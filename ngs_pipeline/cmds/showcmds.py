
from ngs_pipeline import cout, cerr, arg_parser
from ngs_pipeline.cmds import list_commands


def init_argparser():

    p = arg_parser('show all availabel commands')
    return p


def main(args):

    cout("Available commands:")
    for cmd in list_commands():
        cerr(f'  {cmd}')


# EOF
