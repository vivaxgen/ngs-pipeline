#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

__copyright__ = '''
ngs-pl - ngs-pipeline command line executor
[https://github.com/vivaxgen/ngs-pipeline]

(c) 2023-2024 Hidayat Trimarsanto <trimarsanto@gmail.com>

All right reserved.
This software is licensed under MIT license.
Please read the README.txt of this software.
'''

import sys
import os
import logging

# check that we have NGS_PIPELINE_BASE environemt
if 'NGS_PIPELINE_BASE' not in os.environ:
    print('ERROR: please set proper shell enviroment by sourcing activate.sh',
          file=sys.stderr)
    sys.exit(1)

# prepare logging ahead of everything else
if (LOGLEVEL := int(os.environ.get('NGS_PIPELINE_LOGLEVEL', 0))) > 0:
    if (LOGFILE := os.environ.get('NGS_PIPELINE_LOGFILE', '')):
        logging.basicConfig(filename=LOGFILE, level=LOGLEVEL)
    else:
        logging.basicConfig(level=LOGLEVEL)


#from ngs_pipeline import cerr, cexit, subcommands
from ngs_pipeline.subcommands import (_cerr as cerr,
                                      _cexit as cexit,
                                      SubCommands)
import platform


def greet():
    cerr(f'{sys.argv[0].split("/")[-1]} - ngs-pipeline command line interface\n'
         f'[https://github.com/vivaxgen/ngs-pipeline]')
    cerr(f'Host: {platform.uname().node}')


def usage():
    cexit('  usage:\n'
          '      ngs-pl CMD [ARGS]\n'
          '  try: ngs-pl showcmds')


def main():

    cmds = SubCommands(
        modules=['ngs_pipeline.cmds'],
        module_env='NGS_PIPELINE_CMD_MODS',
        env_takes_precedence=True,
        allow_any_script=True,
        allow_shell=True,
        greet_func=greet,
        usage_func=usage,
    )

    cmds.main()


if __name__ == '__main__':
    main()

# EOF
