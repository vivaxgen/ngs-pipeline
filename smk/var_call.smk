
include: "trim_reads.smk"

rule all:
	input:
		"reads/trimmed-{idx}_R1.fastq.gz",
        "reads/trimmed_{idx}_R2.fastq.gz"

# EOF
