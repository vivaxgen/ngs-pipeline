# this is default snakefile for run-targeted-variant-caller

from ngs_pipeline import cexit, cerr
from ngs_pipeline.rules import inc

wf_file = config.get("targeted_variant_caller_wf", "")
if not wf_file:
    cexit(
        "ERROR: please provide the snakefile with"
        " targeted_variant_caller_wf key in config file"
    )

cerr(f"targeted_varcall.smk calls {wf_file}")
include: inc(wf_file)


# EOF