# initialize.smk - ngs-pipeline rules
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

# prepares files for usage

# include utilites.smk and general_params.smk from vivaxGEN ngs-pipeline
include: "utilities.smk"
include: "general_params.smk"


rule wgs:
    input:
        refseq,
        f"{refseq}.fai",
        f"{refseq}.{idx_extension}",
        f"{refseq.removesuffix('.fasta')}.dict",
        strtable_file if strtable_file else [],


rule panelseq:
    input:
        refseq,
        f"{refseq}.fai",
        refmap,


rule snpEff_db:
    input:
        f"{snpEff_data_dir}/{snpEff_db}/snpEffectPredictor.bin",


rule build_snpEff_db:
    input:
        gff_file = f"{snpEff_data_dir}/{snpEff_db}/genes.gff",
        ref_file = f"{snpEff_data_dir}/{snpEff_db}/sequences.fa"
    output:
        snpEff_file = f"{snpEff_data_dir}/{snpEff_db}/snpEffectPredictor.bin",
    shell:
        "snpEff build -c {snpEff_config_file} -dataDir {snpEff_data_dir} -noCheckCds -noCheckProtein -gff3 {snpEff_db}"       


# EOF
