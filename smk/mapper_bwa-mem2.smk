
# BWA-MEM2 mapper

rule reads_mapping:
    input:
        read1 = "trimmed-reads/trimmed-{idx}_R1.fastq.gz",
        read2 = "trimmed-reads/trimmed-{idx}_R2.fastq.gz"
    output:
        bam = "maps/mapped-{idx}.bam",
        stat = "maps/unique_pairs-{idx}.txt.gz"
    log:
        log1 = "logs/bwa-mem2-{idx}.log",
        log2 = "logs/unique_pairs-{idx}.log"
    shell:
        "bwa-mem2 mem -t 16 --secondary=no {refseq} {input.read1} {input.read2} 2> {log.log1}"
        " | samtools fixmate -r -m - {output.bam}"

# EOF
