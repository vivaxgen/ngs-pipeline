# consolidate_samples.py - ngs-pipeline command
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

# to improve the responsiveness during bash autocomplete, do not import heavy
# modules (such as numpy, pandas, etc) here, but instead import them within the
# functions that require the respective heavy modules

import os
import pathlib
import shutil
from ngs_pipeline import (cerr, cexit, arg_parser)

def init_argparser():
    p = arg_parser('move sample directories with uncompleted or unknown state '
                   'to other directory')

    p.add_argument('-o', '--outdir', required=True,
                   help='Output directory for completed samples')
    p.add_argument('-f', '--failed', required=True,
                   help='Target directory for failed samples')
    
    p.add_argument('indirs', nargs='+',
                   help='source directories containing sample directories')
    
    return p


def scan_directory_and_link(a_directory: str | pathlib.Path,
                            completed_dir: str | pathlib.Path,
                            failed_dir: str | pathlib.Path):
    """ return (completed: int, uncompleted: [], unknown: [])"""

    from ngs_pipeline.fileutils import create_relative_symlink

    if not isinstance(a_directory, pathlib.Path):
        a_directory = pathlib.Path(a_directory)

    completed = []
    uncompleted = []
    unknown = []
    for sample_dir in a_directory.iterdir():
        sample_path = pathlib.Path(sample_dir)
        if not sample_path.is_dir():
            continue
        name = sample_path.name
        if (sample_path / '.completed').is_file():
            completed.append((name, sample_path.as_posix()))
            create_relative_symlink(completed_dir / name, sample_path)
        elif (sample_path / '.uncompleted').is_file():
            uncompleted.append((name, sample_dir.as_posix()))
            create_relative_symlink(failed_dir / name, sample_path)
        else:
            unknown.append((name, sample_dir.as_posix()))
            create_relative_symlink(failed_dir / name, sample_path)

    return (completed, uncompleted, unknown)


def consolidate_samples(args):

    completed = []
    uncompleted = []
    unknown = []

    # prepare target directories
    (outdir := pathlib.Path(args.outdir)).mkdir(exist_ok=True)
    (failed := pathlib.Path(args.failed)).mkdir(exist_ok=True)

    for a_dir in args.indirs:
        cerr(f'[Checking sample directory: {a_dir}]')
        compl, uncompl, unkn = scan_directory_and_link(
            a_dir, outdir, failed,
        )
        completed += compl
        uncompleted += uncompl
        unknown += unkn

    cerr('\n===================================='
         ' REPORT '
         '====================================\n')
    if any(uncompleted):
        uncompleted_report = [f'{name}: {path}' for name, path in uncompleted]
        cerr('\n  - '.join(['Uncompleted:'] + uncompleted_report))
        cerr('')
    if any(unknown):
        unknown_report = [f'{name}: {path}' for name, path in unknown]
        cerr('\n  - '.join(['Unknown:'] + unknown_report))
        cerr('')

    cerr(f'[Completed: {len(completed)}, Uncompleted: {len(uncompleted)}, Unknown: {len(unknown)}]')
    cerr('\n')


def main(args):
    consolidate_samples(args)

# EOF
