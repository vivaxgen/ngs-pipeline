
# trimmer_null.smk
#
# use this module if:
# - variant caling for population genetics analysis will be performed (instead of clinical purpose
#   or denovo assembling)
# - the mapper is minimap2 since these mappers will perform soft-clipping ono adapters
# - the variant caller is clair3
# - the model use mv table

from ngs_pipeline.rules import inc

#config['instrument'] = ''
config['correction'] = False
config['libprep'] = 'null'
minlen = 0
maxlen = 0
min_read_quality = 0
min_read_qual = 0
min_avg_quality = 13
min_avg_qual = 13
config["fastplong_cut_tail_window_size"] = 0
config["fastplong_cut_tail_mean_quality"] = 0
config["fastplong_trim_front"] = 0
config["fastplong_trim_tail"] = 0

include: inc("ngs_pipeline::trimmer/fastplong.smk")

# EOF
