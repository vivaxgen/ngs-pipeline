# SPDX-FileCopyrightText: 2026 Hidayat Trimarsanto <trimarsanto@gmail.com>
# SPDX-License-Identifier: MIT

__copyright__ = "(c) 2026 Hidayat (Anto) Trimarsanto <trimarsanto@gmail.com>"
__author__ = "trimarsanto@gmail.com"
__license__ = "MIT"

# this snakemake rules is intended to generate null bam files to allow
# DAG to generate the necessary input files for trimming rules


rule mapping_pe:
    localrule: True
    input:
        read1 = "<sp>trimmed-reads/trimmed-{idx}_R1.fastq.gz",
        read2 = "<sp>trimmed-reads/trimmed-{idx}_R2.fastq.gz",
    output:
        bam = "<sp>maps/mapped-{idx}.bam" if keep_paired_bam else temp("<sp>maps/mapped-{idx}.bam"),
    params:
        sample = get_sample,
    threads: 1
    shell:
        """
        echo -e "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000" | samtools view -b - > {output.bam}
        """




# EOF
