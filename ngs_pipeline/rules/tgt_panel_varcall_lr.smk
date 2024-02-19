# tgt_panel_varcall_lr.smk - ngs-pipeline rules
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2023-2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

from ngs_pipeline import cerr

cerr('Running: tgt_panel_varcall_lr.smk')

include: "utilities.smk"
include: "tgt_params.smk"

if read_files.mode != read_files.ReadMode.SINGLETON:
    raise ValueError('Input files are not single read file per sample(s)')

# prepare sample directory structure
include: "tgt_prepare_sample_files.smk"

# use long-read trimmer:
include: "tgt_trimmer_lr.smk"

# use minimap2 mapper & map merger
include: "tgt_mapper_minimap2_lr.smk"
include: "tgt_merge_map.smk"

# use FreeBayes panel variant-calling
include: "tgt_panel_freebayes.smk"

# include report-generating stuff and report merging
include: "tgt_panel_genreport.smk"


def get_output_file(w):
    print(outdir)
    #return [f'{outdir}/{sample}/vcfs/variants.vcf.gz.csi' for sample in read_files.keys()]
    return [f'{outdir}/{sample}/vcfs/variants.vcf.gz' for sample in read_files.samples()]


rule all:
    input:
        get_output_file

# EOF