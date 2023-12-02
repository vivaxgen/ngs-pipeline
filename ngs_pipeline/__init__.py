
# library providing file common functions

import os
import sys
import argparse
import platform
import argcomplete
import pathlib


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
    return os.environ["NGSENV_BASEDIR"]


def check_NGS_PIPELINE_BASE():
    if 'NGS_PIPELINE_BASE' not in os.environ:
         cexit('ERROR: NGS_PIPELINE_BASE environment is not set. '
              'Please set proper shell enviroment by sourcing relevant activate.sh')
    return os.environ["NGS_PIPELINE_BASE"]     


def get_command_modules():
    """ return a list of modules containing commands"""
    modules = ['ngs_pipeline.cmds']
    if 'NGS_PIPELINE_CMD_MODS' in os.environ:
        additional_modules = [m for m in os.environ['NGS_PIPELINE_CMD_MODS'].split(':') if m]
        return modules + additional_modules
    return modules


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


def get_file_path(filepath: str):
    """ return the actual filepath"""
    if filepath.startswith('/') or filepath.startswith('./') or filepath.startswith('../'):
        return filepath
    return pathlib.Path(check_NGSENV_BASEDIR()) / filepath


def get_snakefile_path(filepath: str, snakefile_root: pathlib.Path):
    """ return the actual snakefile """
    if filepath.startswith('/') or filepath.startswith('./') or filepath.startswith('../'):
        return filepath
    return snakefile_root / filepath


def setup_config(d = {}):
    d['NGSENV_BASEDIR'] = check_NGSENV_BASEDIR()
    d['NGS_PIPELINE_BASE'] = check_NGS_PIPELINE_BASE()
    return d


# EOF
