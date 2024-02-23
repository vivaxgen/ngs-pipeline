# msf_panel_varcall_lr.smk - ngs-pipeline rules
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2023-2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

from ngs_pipeline import cerr

cerr('Running: msf_panel_varcall_lr.smk')

include: "utilities.smk"
include: "msf_params.smk"

if read_files.mode != read_files.ReadMode.SINGLETON:
    raise ValueError('Input files are not single read file per sample(s)')

# prepare sample directory structure
include: "msf_prepare_sample_files.smk"

# use long-read trimmer:
include: "msf_trimmer_lr.smk"

# use minimap2 mapper & map merger
include: "msf_mapper_minimap2_lr.smk"
include: "msf_merge_map.smk"

# use FreeBayes panel variant-calling
include: "msf_varcall_freebayes.smk"

# include report-generating stuff and report merging
include: "msf_panel_genreport.smk"

# include merging vcf
include: "msf_merge_vcf.smk"

def get_individual_output_file(w):
    return [f'{outdir}/{sample}/vcfs/variants.vcf.gz' for sample in read_files.samples()]


cerr(f'Output directory: {outdir}')


rule all:
    input:
        f'{outdir}/merged.vcf.gz'

# EOF
