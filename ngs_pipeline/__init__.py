
# library providing file common functions

import os
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


def check_NGSENV_BASEDIR():
    if 'NGSENV_BASEDIR' not in os.environ:
        cexit('ERROR: NGSENV_BASEDIR environment is not set. '
              'Please set proper shell enviroment by sourcing relevant activate.sh')


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


def add_pgline(alignment_file, pg_dict):
    """ append new PG line using a PG dict, return the full header """

    header = alignment_file.header.to_dict()
    pg_header = header['PG']
    pg_dict['PP'] = pg_header[-1]['ID']
    pg_header.append(pg_dict)
    return header


def get_mode(filename, mode):
    """ return either (r, rb, w, wb) depending on file name """
    if filename.endswith('.bam'):
        return mode + 'b'
    return mode


# EOF
