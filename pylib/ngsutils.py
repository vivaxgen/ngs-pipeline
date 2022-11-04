
# library providing file common functions

import sys
import argparse
import platform
import argcomplete


def cout(msg):
    print(msg, file=sys.stdout)


def cerr(msg):
    print(msg, file=sys.stderr)
    sys.stderr.flush()


def cexit(msg):
    cerr(msg)
    sys.exit(1)


def greet():
    cerr(f'{sys.argv[0].split("/")[-1]} - ngs-pipeline command line interface\n'
         f'[https://github.com/vivaxgen/ngs-pipeline]')
    cerr(f'Host: {platform.uname().node}')


def arg_parser(desc=''):
    p = argparse.ArgumentParser(description=desc)
    p.add_argument('--debug', action='store_true', default=False,
                   help='open ipdb when uncatched exception is raised')
    return p


def run_main(argument_parser=None, main_function=None):

    p = argument_parser() if argument_parser else arg_parser()

    # perform bash autocompletion if needed
    argcomplete.autocomplete(p)

    # parse arguments
    args = p.parse_args()

    # provide greet
    greet()

    if args.debug:
        # run main function under ipdb context
        from ipdb import launch_ipdb_on_exception
        cerr('<< WARN: running in debug mode >>')
        with launch_ipdb_on_exception():
            main_function(args)
    else:
        main_function(args)

# EOF
