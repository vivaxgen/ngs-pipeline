import ngs_pipeline
from ngs_pipeline import cerr, cexit
from ngs_pipeline.snakeutils import (
    get_snakefile_path,
    path_to_str,
    set_default_rule_path,
)

# set the default rule path to ngs_pipeline
set_default_rule_path(ngs_pipeline)


def path(snakefile: str):
    return __path__[0] + "/" + snakefile


__included_snakefiles__ = set()
__void_snakefile__ = path_to_str(
    get_snakefile_path("ngs_pipeline::helper/void.smk", strict_mode=True)
)


def pkg(fn):
    if not fn:
        cexit("include: empty filename")

    fullpath = path_to_str(get_snakefile_path(fn, strict_mode=False))

    if fullpath in __included_snakefiles__:
        cerr(f"including: {fn} (already included, skipping)")
        return __void_snakefile__

    __included_snakefiles__.add(fullpath)
    cerr(f"including: {fn}")
    # cerr(f"fullpath: {fullpath}")
    return fullpath


# alias
inc = pkg

# EOF
