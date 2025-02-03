# msf_mapper_minimap2_pe.smk - ngs-pipeline rules
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2023, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

# required variables:
# - refmap

rule msf_mapping:
    threads: 8
    input:
        read1 = "{pfx}/{sample}/trimmed-reads/trimmed-{idx}_R1.fastq.gz",
        read2 = "{pfx}/{sample}/trimmed-reads/trimmed-{idx}_R2.fastq.gz",
    output:
        bam = "{pfx}/{sample}/maps/{sample}-{idx}.bam",
    log:
        log1 = "{pfx}/{sample}/logs/minimap2-{idx}.log",
        log2 = "{pfx}/{sample}/logs/filter-reads-{idx}.json",
        log3 = "{pfx}/{sample}/logs/filter_reads_region-{idx}.log",
        log4 = "{pfx}/{sample}/logs/fixmate-{idx}.log",
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



# EOF
