
# trimmer_null.smk
#
# use this module if:
# - variant caling for population genetics analysis will be performed (instead of clinical purpose
#   or denovo assembling)
# - the mapper is bwa-mem, bwa-mem2 or minimap since these mappers will perform soft-clipping
#   on the adapters
# - the variant caller is GATK or FreeBayes since these variant-callers will look up at base quality

config['instrument'] = ''
config['correction'] = False
config['libprep'] = 'null'
minlen = 0
maxlen = 0
min_read_quality = 0
min_avg_quality = 13

include: "trimmer_fastp.smk"

