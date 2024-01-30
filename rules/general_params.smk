
# generic parameters

# get base directories
ngs_pipeline_basedir = config['NGS_PIPELINE_BASE']
ngsenv_basedir = config['NGSENV_BASEDIR']

# basic parameters to do processing
refseq = ngsenv_basedir + '/' + config['refseq_file']
ploidy = int(config.get('ploidy', 2))
min_depth = config.get('min_depth', 5)

# thread allocations
thread_allocations = config.get('thread_allocations', {})

# wildcard constraints:
wildcard_constraints:
    sample = r'[.\w-]+',
    idx = '\d',

# EOF
