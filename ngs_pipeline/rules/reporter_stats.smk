# ssf_stats.smk - ngs-pipeline rules
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2024-2026 Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"


rule map_stats:
    threads: 1
    input:
        "<sp>maps/{filename}.bam"
    output:
        "<sp>logs/{filename}.stats.txt"
    params:
        sample = sample,
    shell:
        "samtools stats {input} > {output}"


rule map_samtools_stats:
    # this rule use samtools stats
    threads: 1
    input:
        "<sp>maps/{filename}.bam"
    output:
        "<sp>logs/{filename}.stats.txt-unused"
    params:
        sample = sample,
    shell:
        "samtools stats {input} > {output}"


rule map_sambamba_depth_base:
    # this rule use sambamba depth base
    threads: 1
    input:
        bam = "<sp>maps/{filename}.bam",
        idx = "<sp>maps/{filename}.bam.bai",
    output:
        tsv = "<sp>logs/{filename}.depth-base.tsv.gz",
    log:
        log1 = "<sp>logs/sambamba-depth-{filename}.log",
    params:
        sample = sample,
    shell:
        "sambamba depth base {input.bam} 2> {log.log1} | gzip > {output.tsv}"


rule map_depth:
    # this rule creates stats and depths of any .bam file using samtools stats
    threads: 1
    input:
        '<sp>maps/{filename}.bam'
    output:
        '<sp>logs/{filename}.depths.txt.gz'
    params:
        sample = sample,
    shell:
        'samtools depth {input} | gzip > {output}'


rule depth_plot:
    threads: 1
    input:
        "<sp>logs/mapped-final.depths.txt.gz"
    params:
        sample = sample,
        chroms = ('--chrom ' + ','.join(REGIONS)) if any(REGIONS) else '',
    output:
        '<sp>logs/depths.png'
    log:
        '<sp>logs/plot-depth.txt'
    shell:
        'ngs-pl plot-depth --outplot {output} {params.chroms} --sort '
        '--infile {input} {sample} 2> {log}'


# EOF
