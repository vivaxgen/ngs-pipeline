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
        f"{refseq}",
        f"{refseq}.fai",
        f"{refseq}.bwt.2bit.64",
        f"{refseq.removesuffix('.fasta')}.dict"


rule panelseq:
    input:
        f"{refseq}",
        f"{refseq}.fai",
        f"{refmap}",


# EOF
