
from ngs_pipeline import cerr, fileutils

include: "general_params.smk"

# reference-related configuration
refmap = ngsenv_basedir + '/' + config['refmap_file']
target_regions = ngsenv_basedir + '/' + config['target_regions']
target_variants = ngsenv_basedir + '/' + config['target_variants']
variant_info = ngsenv_basedir + '/' + config.get('variant_info', 'NOFILE')

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
