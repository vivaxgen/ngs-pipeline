
# parameters needed:
# keep_paired_bam = True/False

# BWA / BWA-MEM2 mapper
#
# RG: read groups needs to be defined for proper GATK operation
# ID: unique ID of a collection of reads sequenced together
# SM: sample name
# LB: library prep ID, can be batch or set name for simplicity
# PL: sequencing platform, eg: ILLUMINA

# Anto: not sure why I hardcoded ngs_platform here
ngs_platform="generic"

rule reads_mapping:
    # this rule provides a name-sorted bam file with only mapped paired reads with:
    # - regions mentioned in contaminant_regions removed out, or
    # - regions mentioned in regions filtered-in
    # no further filtering is done (eg. secondary/trans/CPP)
    # the final bam file is suitable for uploading to SRA public database
    threads: thread_allocations.get('bwa_mapping', 4)
    input:
        read1 = "{pfx}/{sample}/trimmed-reads/trimmed-{idx}_R1.fastq.gz",
        read2 = "{pfx}/{sample}/trimmed-reads/trimmed-{idx}_R2.fastq.gz",
        # the following is for sanity check purposes
        refseq = refseq,
        refmap = f"{refseq}.{idx_extension}"

    output:
        #bam = "{pfx}/{sample}/maps/{sample}-{idx}.bam" if keep_paired_bam else temp("{pfx}/maps/{sample}-{idx}.bam"),
        bam = "{pfx}/{sample}/maps/{sample}-{idx}.bam"
    log:
        log0 = "{pfx}/{sample}/logs/reads_mapping-{sample}-{idx}.log",
        log1 = "{pfx}/{sample}/logs/bwa-mem2-{idx}.log",
        log2 = "{pfx}/{sample}/logs/filter-reads-{idx}.json",
        log3 = "{pfx}/{sample}/logs/filter_reads_region-{idx}.log",
        log4 = "{pfx}/{sample}/logs/fixmate-{idx}.log"

    params:
        rg = lambda w: f"-R '@RG\tID:{w.sample}-{w.idx}\tSM:{w.sample}\tLB:LIB-{w.sample}-{w.idx}\tPL:{ngs_platform}'",
        regions = ' '.join(CONTAMINANT_REGIONS) if CONTAMINANT_REGIONS else ' '.join(REGIONS),
        mode = '--remove' if CONTAMINANT_REGIONS else '',
        flags = config.get('bwamem2_flags', ''),
        extra_flags = config.get('bwamem2_extra_flags', ''),
    shell:
        "bwa-mem2 mem -M -t {threads} {params.flags} {params.extra_flags} {params.rg} {refseq} {input.read1} {input.read2} 2> {log.log1}"
        " | ngs-pl filter-reads-region --outstat {log.log2} {params.mode} {params.regions} 2> {log.log3}"
        " | samtools fixmate -m - {output.bam} 2> {log.log4}"

# EOF
