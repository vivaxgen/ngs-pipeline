
import os
from ngs_pipeline import (cerr, cexit, arg_parser,
                          check_NGSENV_BASEDIR, check_NGS_PIPELINE_BASE,
                          get_snakefile_path, setup_config, snakeutils)


def init_argparser():
    p = snakeutils.init_argparser(desc='run arbitrary snakefile')

    # input/output options
    p.add_argument('-o', '--outdir', default='analysis',
                   help='directory for output [analysis/]')
    p.add_argument('infiles', nargs='+',
                   help='FASTQ input files, eg. sample-1.fastq.gz')

    return p


def run_snakefile(args):

    config = dict(infiles=args.infiles, outdir=args.outdir)
    status, elapsed_time = snakeutils.run_snakefile(args, config=config)

    if not status:
        cerr('[WARNING: targeted variant calling did not successfully complete]')
    cerr(f'[Finish targeted variant calling (time: {elapsed_time})]')


def main(args):
    run_snakefile(args)

# EOF
