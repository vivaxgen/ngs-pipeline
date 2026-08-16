# SPDX-FileCopyrightText: 2024-2006 Hidayat (Anto)Trimarsanto <trimarsanto@gmail.com>
# SPDX-License-Identifier: MIT

__copyright__ = "(C) 2024-2006 Hidayat (Anto) Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

# required variables:
# - refmap

rule reads_mapping_lr:
    threads: 8
    input:
        read = "<sp>trimmed-reads/trimmed-{idx}.fastq.gz"
    output:
        #bam = temp("{pfx}/{sample}/maps/{sample}-{idx}.bam")
        bam = temp_unless(get_mapped_bam_file(), keep_paired_bam),
    log:
        log1 = "<sp>logs/minimap2-{idx}.log",
    params:
        rg = lambda w: f"-R @RG\\\\tID:{get_sample(w)}-{w.idx}\\\\tSM:{get_sample(w)}\\\\tLB:LIB-{get_sample(w)}-{w.idx}",
        threads = lambda wildcards, threads: threads - 1,
        flags = config.get('minimap2_flags', ''),
        extra_flags = config.get('minimap2_extra_flags', ''),
    shell:
        "minimap2 -t {params.threads} -a {refmap} -y"
        "  {params.flags} {params.rg} {params.extra_flags}"
        "  {input.read} 2> {log.log1}"
        " | samtools sort -@4 -o {output.bam}"


# EOF

