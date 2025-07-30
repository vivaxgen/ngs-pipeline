# ssf_mapper_bowtie2.smk - ngs-pipeline rule
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2024 Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"


# parameters needed:
# keep_paired_bam = True/False

# bowtie2 is on its own micromamba env, so need to split this into 2 rules:
# the mapping to temp file and the filtering

rule reads_mapping:
    threads: thread_allocations.get('bwa2_mapping', 16)
    input:
        read1 = "trimmed-reads/trimmed-{idx}_R1.fastq.gz",
        read2 = "trimmed-reads/trimmed-{idx}_R2.fastq.gz",
        # the following is for sanity check purposes
        refseq = refseq,
        refmap = f"{refseq}.bt2"

    output:
        bam = temp("maps/temp-{idx}.bam")

    log:
        log0 = f"logs/reads_mapping-{sample}-{{idx}}.log",
        log1 = "logs/bowtie2-{idx}.log",
        log2 = "logs/filter-reads-{idx}.json",
        log3 = "logs/filter_reads_region-{idx}.log",
        log4 = "logs/fixmate-{idx}.log",

    params:
        sample = sample,
        rg = lambda w: f"-R '@RG\tID:{sample}-{w.idx}\tSM:{sample}\tLB:LIB-{sample}-{w.idx}\tPL:{ngs_platform}'",
        regions = ' '.join(CONTAMINANT_REGIONS) if CONTAMINANT_REGIONS else ' '.join(REGIONS),
        mode = '--remove' if CONTAMINANT_REGIONS else '',
        flags = config.get('bowtie2_flags', ''),
        extra_flags = config.get('bowtie2_extra_flags', ''),

    shell:
        "{shell_prefix_bowtie2}"
        "bowtie2  -M -t {threads} {params.flags} {params.extra_flags} {params.rg} {refseq} {input.read1} {input.read2} 2> {log.log1}"
        " | ngs-pl filter-reads-region --outstat {log.log2} {params.mode} {params.regions} 2> {log.log3}"
        " | samtools fixmate -m - {output.bam} 2> {log.log4}"

rule reads_mapping_filtering:
    threads: 2
    input:
        bam = "maps/temp-{idx}.sam"
    output:
        bam = "maps/mapped-{idx}.bam" if keep_paired_bam else temp("maps/mapped-{idx}.bam"),

    shell:
        'cat {input} '
        " | ngs-pl filter-reads-region --outstat {log.log2} {params.mode} {params.regions} 2> {log.log3}"
        " | samtools fixmate -m - {output.bam} 2> {log.log4}"

# EOF
