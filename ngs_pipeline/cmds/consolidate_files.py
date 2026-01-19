# check_configs.py - ngs-pipeline command
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

# to improve the responsiveness during bash autocomplete, do not import heavy
# modules (such as numpy, pandas, etc) here, but instead import them within the
# functions that require the respective heavy modules


import os
import pathlib
import shutil
from ngs_pipeline import arg_parser, cerr


filetypes = {
    "map": ("maps/mapped-final.bam", "{sample}.bam"),
    "depth-base": ("logs/mapped-final.depth-base.tsv.gz", "{sample}.depth-base.tsv.gz"),
}


def init_argparser():
    p = arg_parser()

    p.add_argument("-o", "--outdir", help="output directory")
    p.add_argument(
        "--src", required=True, help="source file path relative to sample directory"
    )
    p.add_argument(
        "--dst",
        required=True,
        help="destination file name format, use {sample} to indicate sample name",
    )
    p.add_argument(
        "--copy", action="store_true", help="copy files instead of creating symlinks"
    )
    p.add_argument(
        "--hardlink", action="store_true", help="create hardlinks instead of symlinks"
    )
    p.add_argument("indirs", nargs="+")
    return p


def consolidate_files(args):

    from ngs_pipeline.fileutils import create_relative_symlink

    src, dst = args.src, args.dst
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    total_dir = 0
    consolidated = 0
    failed = 0
    for indir in args.indirs:
        indir = pathlib.Path(indir)
        for sampledir in indir.iterdir():
            if not sampledir.is_dir():
                cerr("Skipping non-directory: {sampledir}")
                continue
            total_dir += 1
            sample = sampledir.name
            srcfile = sampledir / src
            if not srcfile.exists():
                cerr(f"Source file not found: {srcfile}")
                failed += 1
                continue
            dstfile = outdir / dst.format(sample=sample)
            if args.copy:
                shutil.copyfile(srcfile, dstfile, follow_symlinks=True)
            elif args.hardlink:
                os.link(srcfile, dstfile)
            else:
                create_relative_symlink(dstfile, srcfile, force=True)
            consolidated += 1

    cerr(
        f"[Consolidated files [{src}]: {consolidated}, failed: {failed} "
        f"from total: {total_dir} sample directories]"
    )


def main(args):
    consolidate_files(args)


# EOF
