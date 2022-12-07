
import os
import pathlib

ngs_pipeline_basedir = os.environ['NGS_PIPELINE_BASE']
ngsenv_basedir = os.environ['NGSENV_BASEDIR']

refmap = ngsenv_basedir + '/' + config.get('refmap_file', 'NOFILE')
refseq = ngsenv_basedir + '/' + config['refseq_file']
knownsites_file = ngsenv_basedir + '/' + config.get('knownsites_file', '')

java_opts = ''

# check available read files

platform = 'ILLUMINA' if config.get('instrument', '').lower() in ['miseq', 'nextseq', 'novaseq'] else config['platform']

REGIONS = config.get('regions', [])

wildcard_constraints:
    idx = '\\d+'

# EOF
