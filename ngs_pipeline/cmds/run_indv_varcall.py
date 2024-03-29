__copyright__ = '''
run_indv_varcall.py - ngs-pipeline command line
[https://github.com/vivaxgen/ngs-pipeline]

(c) 2022-2024 Hidayat Trimarsanto <trimarsanto@gmail.com>

All right reserved.
This software is licensed under MIT license.
Please read the README.txt of this software.
'''

# to improve the responsiveness during bash autocomplete, do not import heavy
# modules (such as numpy, pandas, etc) here, but instead import them within the
# functions that requires respective heavy modules

import sys
import os
from ngs_pipeline import cerr, cexit, snakeutils


# usage: run_varcall.py

def init_argparser():
    p = snakeutils.init_argparser(desc='run individual sample variant caller')
    p.arg_dict['snakefile'].choices = ['var_call.smk']
    p.arg_dict['snakefile'].default = 'var_call.smk'
    p.add_argument('target')
    p.add_argument('sample')

    return p


def run_indv_varcall(args):

    import pathlib

    # simple sanity checks
    cwd = pathlib.Path.cwd()
    sample = cwd.parts[-1]
    if not (cwd / 'reads').is_dir():
        cexit(f'ERROR: current directory {cwd} does not have "reads" directory!')

    # remove .completed or .uncompleted screen first if exists
    (cwd / '.completed').unlink(missing_ok=True)
    (cwd / '.uncompleted').unlink(missing_ok=True)

    status, elapsed_time = snakeutils.run_snakefile(args)
    open(cwd / ('.uncompleted' if not status else '.completed'), 'w').write(f'{elapsed_time}')


def main(args):
    run_indv_varcall(args)

# EOF
