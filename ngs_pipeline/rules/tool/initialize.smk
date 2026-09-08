# initialize.smk - ngs-pipeline rules
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"


from ngs_pipeline.rules import pkg
# prepares files for usage


# include utilites.smk and general_params.smk from vivaxGEN ngs-pipeline
include: pkg("ngs_pipeline::helper/utilities.smk")
include: pkg("ngs_pipeline::general_params.smk")


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

all_variant_vcf = [config.get(k) for k in config.keys() if k.startswith("target_variants_vcf")]

rule variant_vcf:
    input:
        *[f"{get_abspath(vcf)}.csi" for vcf in all_variant_vcf]

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
