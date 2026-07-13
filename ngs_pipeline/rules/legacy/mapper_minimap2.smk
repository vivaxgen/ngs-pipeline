
# mapper_minimap2_pe.smk - ngs-pipeline rules
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2023-2025, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

# this is unified msf (multi-sample flow) and ssf (single-sample flow)  utilizing  <path-var> and <sample_dir> for msf.
# The rule is designed to be flexible and can be used in both flows without modification. The required variables are
# defined in the config file and can be accessed using the config object.

# required variables:
# - refmap

rule ssf_mapping:
    threads: 8
    input:
        read1 = "trimmed-reads/trimmed-{idx}_R1.fastq.gz",
        read2 = "trimmed-reads/trimmed-{idx}_R2.fastq.gz",
    output:
        bam = f"maps/{sample}-{{idx}}.bam",
    log:
        log1 = "logs/minimap2-{idx}.log",
        log2 = "logs/filter-reads-{idx}.json",
        log3 = "logs/filter_reads_region-{idx}.log",
        log4 = "logs/fixmate-{idx}.log",


rule msf_mapping:
    threads: 8
    input:
        read1 = "<sample_dir>/trimmed-reads/trimmed-{idx}_R1.fastq.gz",
        read2 = "<sample_dir/trimmed-reads/trimmed-{idx}_R2.fastq.gz",
        meta = "<sample_dir>/{sample}" 
    output:
        bam = "<sample_dir>/maps/{sample}-{idx}.bam",
    log:
        log1 = "<sample_dir>/logs/minimap2-{idx}.log",
        log2 = "<sample_dir>/logs/filter-reads-{idx}.json",
        log3 = "<sample_dir>/logs/filter_reads_region-{idx}.log",
        log4 = "<sample_dir>/logs/fixmate-{idx}.log",
    params:
        rg = lambda w: f"-R @RG\\\\tID:{w.sample}-{w.idx}\\\\tSM:{w.sample}\\\\tLB:LIB-{w.sample}-{w.idx}",
        threads = lambda wildcards, threads: threads - 1,
        regions = ' '.join(CONTAMINANT_REGIONS) if CONTAMINANT_REGIONS else ' '.join(REGIONS),
        mode = '--remove' if CONTAMINANT_REGIONS else '',
        flags = config.get('minimap2_flags', ''),
        extra_flags = config.get('minimap2_extra_flags', ''),
    shell:
        "minimap2 -t {params.threads} -a {refmap} {params.rg}"
        "  {params.flags} {params.extra_flags}"
        "  {input.read1} {input.read2} 2> {log.log1}"
        " | samtools collate -u -O -"
        " | samtools fixmate -m - - 2> {log.log4}"
        " | ngs-pl filter-reads-region -o {output.bam} --outstat {log.log2} {params.mode} {params.regions} 2> {log.log3}"



rule reads_mapping:
    threads: 8
    input:
        read1 = "trimmed-reads/trimmed-{idx}_R1.fastq.gz",
        read2 = "trimmed-reads/trimmed-{idx}_R2.fastq.gz"
    output:
        bam = "maps/mapped-{idx}.bam",
        stat = "maps/unique_pairs-{idx}.txt.gz"
    log:
        log1 = "logs/minimap2-{idx}.log",
        log2 = "logs/unique_pairs-{idx}.log"
    shell:
        "minimap2 -ax sr -t {threads} --secondary=no {refmap} {input.read1} {input.read2} 2> {log.log1}"
        " | {ngs_pipeline_basedir}/bin/filter_uniquepair.py --min_match_len 23 --max_nm 0.29 --outstat {output.stat} -o - - 2> {log.log2}"
        " | samtools fixmate -r -m - {output.bam}"

# EOF
