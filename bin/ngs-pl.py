#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

__copyright__ = '''
ngs-pl - ngs-pipeline command line executor
[https://github.com/vivaxgen/ngs-pipeline]

(c) 2023 Hidayat Trimarsanto <trimarsanto@gmail.com>

All right reserved.
This software is licensed under MIT license.
Please read the README.txt of this software.
'''

import sys
import os

# check that we have NGS_PIPELINE_BASE environemt
if 'NGS_PIPELINE_BASE' not in os.environ:
    print('ERROR: please set proper shell enviroment by sourcing activate.sh',
          file=sys.stderr)
    sys.exit(1)


from ngs_pipeline import cerr, cexit, run_main
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

    tokens = []

    if '_ARGCOMPLETE' in os.environ:
        # we are invoked from argcomplete
        line = os.environ.get('COMP_LINE', '')
        tokens = line.split()
        if len(tokens) == 1 or (len(tokens) == 2 and not line.endswith(' ')):
            autocomplete(tokens)
        os.environ['COMP_LINE'] = line.split(' ', 1)[1]
        os.environ['COMP_POINT'] = str(len(os.environ['COMP_LINE']))
        cmd = tokens[1]
    else:
        greet()
        if len(sys.argv) == 1:
            usage()
        cmd = sys.argv[1]

    from ngs_pipeline.cmds import run_main
    run_main(tokens[1:] if any(tokens) else sys.argv[1:])


def autocomplete(tokens):

    from ngs_pipeline.cmds import list_commands

    # prepare line

    last_token = tokens[-1]

    # prepare the completion lists
    cmds = list_commands()

    if len(tokens) == 1:
        completions = sorted(cmds)

    else:
        completions = sorted(
            [opt for opt in cmds if opt.startswith(last_token)]
        )

    # send results through fd 8
    ifs = os.environ.get('IFS', '\013')
    out_stream = os.fdopen(8, 'w')
    out_stream.write(ifs.join(completions))
    sys.exit(0)


if __name__ == '__main__':
    main()

# EOF
