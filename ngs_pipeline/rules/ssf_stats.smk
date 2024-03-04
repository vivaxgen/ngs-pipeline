# ssf_stats.smk - ngs-pipeline rules
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2024 Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"


rule map_samtools_stats:
    # this rule use samtools stats
    threads: 1
    input:
        "maps/{filename}.bam"
    output:
        "logs/{filename}.stats.txt-unused"
    params:
        sample = sample,
    shell:
        "samtools stats {input} > {output}"


rule map_sambamba_depth_base:
    # this rule use sambamba depth base
    threads: 1
    input:
        "maps/{filename}.bam"
    output:
        "logs/{filename}.depth-base.tsv.gz"
    params:
        sample = sample,
    shell:
        "sambamba depth base {input} | gzip > {output}"

# EOF
