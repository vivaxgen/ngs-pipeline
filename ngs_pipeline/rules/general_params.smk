
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
    idx = r'\d',

# for indexing, use the correct extension for bwa / bwa-mem2

bwa_bin = config.get('bwa_bin', 'bwa-mem2')
if bwa_bin == 'bwa-mem2':
    idx_extension = 'bwt.2bit.64'
elif bwa_bin == 'bwa':
    idx_extension = 'bwt'
else:
    raise RuntimeError(f'ERR: the bwa_bin {bwa_bin} is not recognized')

# EOF
