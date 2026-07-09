from ngs_pipeline import cerr


def path(snakefile: str):
    return __path__[0] + "/" + snakefile


def inc(fn):
    cerr(f"including: {fn}")
    return fn


# EOF
