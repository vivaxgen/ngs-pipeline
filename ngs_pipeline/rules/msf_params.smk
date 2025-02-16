# msf_params.smk - ngs-pipeline rules
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2023, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

from ngs_pipeline import cerr, fileutils

include: "general_params.smk"
include: "params_region.smk"

# reference-related configuration
target_regions = get_abspath(config['target_regions'], ngsenv_basedir)
target_variants = get_abspath(config['target_variants'], ngsenv_basedir)
variant_info = get_abspath(config.get('variant_info', 'NOFILE'), ngsenv_basedir)

# input/output related configuration
read_files = fileutils.ReadFileDict(config['infiles'], config['underscore'])
outdir = config['outdir']

# read quality-related configuration
min_read_qual = config['min_read_qual']
max_read_len = config['max_read_len']
min_read_len = config['min_read_len']
headcrop = config.get('headcrop', 0)
tailcrop = config.get('tailcrop', 0)

# variant quality-related coniguration
min_variant_qual = config.get('min_variant_qual', 30)

# required for some snakemake slurm profile for naming the slurm job
sample = 'all'

# EOF
