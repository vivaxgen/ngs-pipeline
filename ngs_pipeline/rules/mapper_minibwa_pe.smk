# SPDX-FileCopyrightText: 2026 Hidayat Trimarsanto <trimarsanto@gmail.com>
# SPDX-License-Identifier: MIT

__copyright__ = "(c) 2026 Hidayat (Anto) Trimarsanto <trimarsanto@gmail.com>"
__author__ = "trimarsanto@gmail.com"
__license__ = "MIT"


rule reads_mapping:
    threads: 8
    input:
        read1 = "<sp>trimmed-reads/trimmed-{idx}_R1.fastq.gz",
        read2 = "<sp>trimmed-reads/trimmed-{idx}_R2.fastq.gz",
    output:
        bam = temp_unless(get_mapped_bam_file(), keep_paired_bam),
    log:
        log1 = "<sp>logs/minibwa-{idx}.log",
        log2 = "<sp>logs/filter-reads-{idx}.json",
        log3 = "<sp>logs/filter_reads_region-{idx}.log",
        log4 = "<sp>logs/fixmate-{idx}.log",
    params:
        sample = get_sample,
        rg = lambda w: f"-R '@RG\\tID:{get_sample(w)}-{w.idx}\\tSM:{get_sample(w)}\\tLB:LIB-{get_sample(w)}-{w.idx}\\tPL:{ngs_platform}'",
        threads = lambda wildcards, threads: threads - 1,
        regions = ' '.join(CONTAMINANT_REGIONS) if CONTAMINANT_REGIONS else ' '.join(REGIONS),
        mode = '--remove' if CONTAMINANT_REGIONS else '',
        flags = config.get('minibwa_flags', ''),
        extra_flags = config.get('minibwa_extra_flags', ''),
    shell:
        "minibwa map -t {params.threads} {params.rg}"
        "  {params.flags} {params.extra_flags}"
        "  {refseq} {input.read1} {input.read2} 2> {log.log1}"
        #" | samtools collate -u -O -"
        " | samtools fixmate -m - - 2> {log.log4}"
        " | ngs-pl filter-reads-region -o {output.bam} --outstat {log.log2} {params.mode} {params.regions} 2> {log.log3}"

