
# BWA-MEM2 mapper
#
# RG: read groups needs to be defined for proper GATK operation
# ID: unique ID of a collection of reads sequenced together
# SM: sample name
# LB: library prep ID, can be batch or set name for simplicity
# PL: sequencing platform, eg: ILLUMINA


rule reads_mapping:
    threads: 16
    input:
        read1 = "trimmed-reads/trimmed-{idx}_R1.fastq.gz",
        read2 = "trimmed-reads/trimmed-{idx}_R2.fastq.gz"
    output:
        bam = temp("maps/mapped-{idx}.bam"),
    log:
        log1 = "logs/bwa-mem2-{idx}.log",
    params:
        rg = lambda w: f"-R '@RG\tID:{sample}-{w.idx}\tSM:{sample}\tLB:LIB-{sample}-{w.idx}\tPL:{platform}'",
        regions = ' '.join(REGIONS)
    shell:
        "bwa-mem2 mem -M -t {threads} {params.rg} {refseq} {input.read1} {input.read2} 2> {log.log1}"
        " | filter_reads.py {params.regions}"
        " | samtools fixmate -r -m - {output.bam}"

# EOF
