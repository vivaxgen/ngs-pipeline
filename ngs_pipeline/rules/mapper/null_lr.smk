# SPDX-FileCopyrightText: 2026 Hidayat Trimarsanto <trimarsanto@gmail.com>
# SPDX-License-Identifier: MIT

__copyright__ = "(c) 2026 Hidayat (Anto) Trimarsanto <trimarsanto@gmail.com>"
__author__ = "trimarsanto@gmail.com"
__license__ = "MIT"

# this snakemake rules is intended to generate null bam files to allow
# DAG to generate the necessary input files for trimming rules

include: "mapper_null_pe.smk"

use rule mapping_pe as mapping_lr:
    input:
        read0 = "<sp>trimmed-reads/trimmed-{idx}.fastq.gz",

# EOF
