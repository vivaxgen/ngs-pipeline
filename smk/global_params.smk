
import os
import pathlib

ngs_pipeline_basedir = os.environ['NGS_PIPELINE_BASE']
ngsenv_basedir = os.environ['NGSENV_BASEDIR']

refmap = ngsenv_basedir + '/' + config.get('refmap_file', 'NOFILE')
refseq = ngsenv_basedir + '/' + config['refseq_file']
knownsites_file = ngsenv_basedir + '/' + config.get('knownsites_file', '')
ploidy = int(config.get('ploidy', 2))

# parameters to save storage space
keep_paired_bam = config.get('keep_paired_bam', False)
keep_proper_bam = config.get('keep_proper_bam', False)
keep_filtered_bam = config.get('keep_filtered_bam', False)
keep_deduplicated_bam = config.get('keep_deduplicated_bam', False)
keep_recalibrated_bam = config.get('keep_recalibrated_bam', False)

# parameters for read trimming
minlen = int(config['minlen']) if 'minlen' in config else int(config['read_length'] / 3)
maxlen = int(config['maxlen']) if 'maxlen' in config else 0
min_quality = int(config.get('min_read_qual', 15))
min_avg_quality = int(config.get('min_avg_qual', 0))
correction = config.get('correction', False)
instrument = config.get('instrument', None)
read_filters = config.get('read_filters', '')


java_opts = ''

# check available read files

platform = 'ILLUMINA' if config.get('instrument', '').lower() in ['miseq', 'nextseq', 'novaseq'] else config['platform']

PARTIALS = None
REGIONS = config.get('regions', [])
if isinstance(REGIONS, dict):
    PARTIALS = REGIONS
    REGIONS = PARTIALS.keys()

CONTAMINANT_REGIONS = config.get('contaminant_regions', [])

wildcard_constraints:
    idx = '\\d+',
    sample = '[\\w-]+'

# EOF
