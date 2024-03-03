# move_failed_samples.py - ngs-pipeline command
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

# to improve the responsiveness during bash autocomplete, do not import heavy
# modules (such as numpy, pandas, etc) here, but instead import them within the
# functions that require the respective heavy modules

import os
import pathlib
import shutil
from ngs_pipeline import (cerr, cexit, arg_parser, get_snakefile_path,
                          snakeutils)

def init_argparser():
    p = arg_parser('move sample directories with uncompleted or unknown state '
                   'to other directory')

    p.add_argument('-o', '--outdir', required=True,
                   help='Out or target diretory to remove to')
    
    p.add_argument('indirs', nargs='+',
                   help='source directories containing sample directories')
    
    return p


def scan_directory_and_move(a_directory: str | pathlib.Path,
                            target_dir: str | pathlib.Path):
    """ return (completed: int, uncompleted: [], unknown: [])"""

    if not isinstance(a_directory, pathlib.Path):
        a_directory = pathlib.Path(a_directory)

    completed = 0
    uncompleted = []
    unknown = []
    for sample_dir in a_directory.iterdir():
        sample_path = pathlib.Path(sample_dir)
        if not sample_path.is_dir():
            continue
        if (sample_path / '.completed').is_file():
            completed += 1
        elif (sample_path / '.uncompleted').is_file():
            uncompleted.append(sample_dir.as_posix())
            shutil.move(sample_path, target_dir)
        else:
            unknown.append(sample_dir.as_posix())
            shutil.move(sample_path, target_dir)

    return (completed, uncompleted, unknown)


def move_failed_samples(args):

    completed = 0
    uncompleted = []
    unknown = []

    for a_dir in args.indirs:
        cerr(f'[Checking sample directory: {a_dir}]')
        compl, uncompl, unkn = scan_directory_and_move(
            a_dir, args.outdir
        )
        completed += compl
        uncompleted += uncompl
        unknown += unkn

    cerr('\n===================================='
         ' REPORT '
         '====================================\n')
    if any(uncompleted):
        cerr('\n  - '.join(['Uncompleted:'] + uncompleted))
        cerr('')
    if any(unknown):
        cerr('\n  - '.join(['Unknown:'] + unknown))
        cerr('')

    cerr(f'[Completed: {completed}, Uncompleted: {len(uncompleted)}, Unknown: {len(unknown)}]')
    cerr('\n')


def main(args):
    move_failed_samples(args)

# EOF
