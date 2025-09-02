# version.py - ngs-pipeline command line
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(c) 2025, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

# to improve the responsiveness during bash autocomplete, do not import heavy
# modules (such as numpy, pandas, etc) here, but instead import them within the
# functions that require the respective heavy modules

from ngs_pipeline import cout, cerr, arg_parser, check_NGS_PIPELINE_BASE


def init_argparser():

    p = arg_parser("show version")
    p.add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help="show git tag using git log command",
    )
    return p


def get_git_hash(repo_dir, label, verbose=False):
    """get git has from a repo_dir"""
    import pathlib
    import subprocess

    git_dir = pathlib.Path(repo_dir) / ".git"
    if not git_dir.exists():
        git_hash = "N/A"
    elif verbose:
        output = subprocess.check_output(
            ["git", "--git-dir", git_dir, "log", "-1", "--source"]
        ).decode()
        git_hash = output.splitlines()[0]
    else:
        git_ref_line = open(git_dir / "HEAD").read().strip().split()
        if len(git_ref_line) != 2:
            return ""
        git_hash = open(git_dir / git_ref_line[1]).read().strip()
    return f"{label} (git hash): {git_hash}"


def version(args):

    git_hash = get_git_hash(
        check_NGS_PIPELINE_BASE(), "NGS-Pipeline", verbose=args.verbose
    )
    cout(git_hash)


def main(args):
    version(args)


# EOF
