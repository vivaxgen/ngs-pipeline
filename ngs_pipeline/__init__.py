# __init__.py - ngs-pipeline
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(c) 2023, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

# library providing file common functions

import os
import sys
import argparse
import platform
import argcomplete
import pathlib
import types


def cout(msg):
    print(msg, file=sys.stdout)


def cerr(msg):
    print(msg, file=sys.stderr)
    sys.stderr.flush()


def cexit(msg, err_code=1):
    cerr(msg)
    sys.exit(err_code)


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


def is_abs_or_rel_path(filepath: str | pathlib.Path):
    filepath = (filepath.as_posix() if isinstance(filepath, pathlib.Path)
                else filepath)
    if (
        filepath.startswith('/')
        or filepath.startswith('./')
        or filepath.startswith('../')
    ):
        return True
    return False


def get_file_path(filepath: str):
    """ return the actual filepath"""
    if is_abs_or_rel_path(filepath):
        return filepath
    return pathlib.Path(check_NGSENV_BASEDIR()) / filepath


def get_snakefile_path(filepath: str | pathlib.Path,
                       snakefile_root: pathlib.Path | None = None,
                       from_module: types.ModuleType | None = None):
    """ return the actual snakefile """
    filepath_str = filepath
    if isinstance(filepath, pathlib.Path):
        filepath_str = filepath.as_posix()
    if is_abs_or_rel_path(filepath):
        return filepath
    if from_module is not None:
        snakefile_root = pathlib.Path(from_module.__path__[0]) / 'rules'
    return snakefile_root / filepath


def setup_config(d = {}):
    d['NGSENV_BASEDIR'] = check_NGSENV_BASEDIR()
    d['NGS_PIPELINE_BASE'] = check_NGS_PIPELINE_BASE()
    return d


# EOF
