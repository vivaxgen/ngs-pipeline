# check_configs.py - ngs-pipeline command
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

# to improve the responsiveness during bash autocomplete, do not import heavy
# modules (such as numpy, pandas, etc) here, but instead import them within the
# functions that require the respective heavy modules


import os
import pathlib
from ngs_pipeline import arg_parser


filetypes = {
    'map': ('maps/mapped-final.bam',
            '{sample}.bam'),
    'depth-base': ('logs/mapped-final.depth-base.tsv.gz',
                   '{sample}.depth-base.tsv.gz'),
}


def init_argparser():
    p = arg_parser()

    p.add_argument('-o', '--outdir',
                   help='output directory')
    p.add_argument('-t', '--type',
                   choices = filetypes.keys(),
                   help='type of file to consolidate')
    p.add_argument('indirs', nargs='+')
    return p


def consolidate_reports(args):

    from ngs_pipeline.fileutils import create_relative_symlink

    src, dst = filetypes[args.type]
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    for indir in args.indirs:
        indir = pathlib.Path(indir)
        for sampledir in indir.iterdir():
            sample = sampledir.name
            srcfile = sampledir / src
            dstfile = outdir / dst.format(sample=sample)
            create_relative_symlink(dstfile, srcfile)


def main(args):
    consolidate_reports(args)

# EOF
