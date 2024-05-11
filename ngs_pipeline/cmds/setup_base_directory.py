# initialize.py - ngs-pipeline command
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

# to improve the responsiveness during bash autocomplete, do not import heavy
# modules (such as numpy, pandas, etc) here, but instead import them within the
# functions that require the respective heavy modules

import os
from pathlib import Path
import ngs_pipeline


def init_argparser():

    p = ngs_pipeline.arg_parser('generate base environment directory')
    p.add_argument('--rcfile', default='')
    p.add_argument('--ngs-pipeline-base', default='')
    p.add_argument('basedir')

    return p


# below is activate script in the NGSENV_BASEDIR directory
# which will be directory-specific
activation_script_template = '''#!/usr/bin/env bash

# -- user-modified variables --
export BASHRC={BASHRC}
export NGS_PIPELINE_BASE={NGS_PIPELINE_BASE}

# -- end of user-modified variables --

. ${{BASHRC}}

# get the directory where this script resides
_script="$(readlink -f ${{BASH_SOURCE[0]}})"
_mydir="$(dirname $_script)"
export NGSENV_BASEDIR=${{_mydir}}

if [[ "${{BASH_SOURCE[0]}}" == "${{0}}" ]]; then
    # new shell will be spawned
    export SPAWN_SHELL=1
fi

. ${{NGS_PIPELINE_BASE}}/bin/activate


'''

# below is the bashrc script in the NGSENV_BASEDIR directory
# that user can further customize
bashrc_source_template = '''
# put additional settings here

# EOF
'''

def setup_base_directory(args):
    
    # create directory structure
    basedir = Path(args.basedir)
    (basedir / 'configs' / 'refs').mkdir(exist_ok=True, parents=True)
    (basedir / 'sets').mkdir(exist_ok=True, parents=True)

    # check if we are under vivaxGEN base install

    VVG_BASEDIR = os.environ.get('VVG_BASEDIR', '')
    if VVG_BASEDIR:
        # we use vivaxGEN vvg-base to generate script file

        import subprocess
        VVGBIN = os.environ.get('VVGBIN')
        subprocess.call([
            f'{VVGBIN}/generate-activation-script.py',
            '-o', (basedir / 'activate').as_posix(),
            '-b', VVG_BASEDIR,
            '-e', f'NGSENV_BASEDIR={basedir.resolve()} ',
            '-e', f'NGS_PROMPT={basedir.name}'
        ])
        
    else:

        if not args.ngs_pipeline_base:
            args.ngs_pipeline_base = ngs_pipeline.check_NGS_PIPELINE_BASE()

        activation_file = basedir / 'activate'

        with open(activation_file, 'w') as f_out:
            f_out.write(
                activation_script_template.format(
                    BASHRC=args.rcfile,
                    NGS_PIPELINE_BASE=args.ngs_pipeline_base
                )
            )
        activation_file.chmod(0o775)

    with open(basedir / 'configs' / 'bashrc', 'w') as f_out:
        f_out.write(bashrc_source_template)


def main(args):
    setup_base_directory(args)


# EOF
