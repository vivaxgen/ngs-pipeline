
rule reads_mapping:
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
        "minimap2 -ax sr -t 16 --secondary=no {refmap} {input.read1} {input.read2} 2> {log.log1}"
        " | {ngs_pipeline_basedir}/bin/filter_uniquepair.py --min_match_len 23 --max_nm 0.29 --outstat {output.stat} -o - - 2> {log.log2}"
        " | samtools fixmate -r -m - {output.bam}"

# EOF
