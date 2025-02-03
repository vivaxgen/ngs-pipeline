# msf_panel_varcall_pe.smk - ngs-pipeline rules
# [https://github.com/vivaxgen/ngs-pipeline]

__author__ = "Hidayat Trimarsanto"
__copyright__ = "(C) 2024, Hidayat Trimarsanto"
__email__ = "trimarsanto@gmail.com,hidayat.trimarsanto@menzies.edu.au"
__license__ = "MIT"

# this snakemake rules is intended for panel variant calling of pair-end reads

from ngs_pipeline import cerr

cerr('Running: msf_panel_varcall_lr.smk')

include: "utilities.smk"
include: "msf_params.smk"

if read_files.mode != read_files.ReadMode.PAIRED_END:
    raise ValueError('Input files are not paired-end read file per sample(s)')

# prepare sample directory structure
include: "msf_prepare_sample_files.smk"

# use null trimmer since we rely on minimap2 to perform soft-clipping
# on primers and adapters:
include: "msf_trimmer_null.smk"

# use minimap2 mapper & map merger
include: "msf_mapper_minimap2_pe.smk"
include: "msf_merge_map.smk"

# use FreeBayes panel variant-calling as default
include: config.get('msf_varcall_wf', 'msf_varcall_freebayes.smk')

# include report-generating stuff and report merging
include: "msf_panel_genreport.smk"

# include merging vcf
include: "msf_merge_vcf.smk"

def get_individual_output_file(w):
    return [f'{outdir}/samples/{sample}/vcfs/variants.vcf.gz' for sample in read_files.samples()]


cerr(f'Output directory: {outdir}')


rule all:
    input:
        f'{outdir}/merged.vcf.gz'


# EOF
