
import os
from ngs_pipeline import (cerr, cexit, arg_parser,
                          check_NGSENV_BASEDIR, check_NGS_PIPELINE_BASE,
                          get_snakefile_path, setup_config, snakeutils)


def init_argparser(desc=None):
    p = snakeutils.init_argparser(desc=desc or 'run arbitrary snakefile')
    return p


def run_snakefile(args, config: dict = {}):

    status, elapsed_time = snakeutils.run_snakefile(
        args,
        config=setup_config(config)
    )

    if not status:
        cerr(f'[WARNING: snakefile {args.snakefile} did not successfully '
             'complete]')
    cerr(f'[Finish running snakefile {args.snakefile} (time: {elapsed_time})]')


def main(args):
    run_snakefile(args)

# EOF
