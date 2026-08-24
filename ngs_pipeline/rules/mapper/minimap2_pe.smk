# msf_mapper_minimap2_pe.smk - ngs-pipeline rules
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2023-2026 Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

# required variables:
# - refmap

include: "init.smk"

rule reads_mapping_pe:
    threads: 8
    input:
        read1 = "{anypath}trimmed-reads/trimmed-{idx}_R1.fastq.gz",
        read2 = "{anypath}trimmed-reads/trimmed-{idx}_R2.fastq.gz",
    output:
        #bam = temp_unless(replace_sp("<sp>maps/{sample}-{idx}.bam"), keep_paired_bam),
        bam = temp_unless("{anypath}maps/mapped-{idx}.bam", keep_paired_bam),
    wildcard_constraints:
        sample = r'[.\w-]+'  # Explicitly overriding locally fixes the parsing bug
    log:
        log1 = "{anypath}logs/minimap2-{idx}.log",
        log2 = "{anypath}logs/filter-reads-{idx}.json",
        log3 = "{anypath}logs/filter_reads_region-{idx}.log",
        log4 = "{anypath}logs/fixmate-{idx}.log",
    params:
        sample = get_sample,
        rg = lambda w: f"-R @RG\\\\tID:{get_sample(w)}-{w.idx}\\\\tSM:{get_sample(w)}\\\\tLB:LIB-{get_sample(w)}-{w.idx}",
        threads = lambda wildcards, threads: threads - 1,
        regions = ' '.join(CONTAMINANT_REGIONS) if CONTAMINANT_REGIONS else ' '.join(REGIONS),
        mode = '--remove' if CONTAMINANT_REGIONS else '',
        flags = config.get('minimap2_flags', ''),
        extra_flags = config.get('minimap2_extra_flags', ''),
    shell:
        "minimap2 -t {params.threads} -ax sr --MD {refmap} {params.rg}"
        "  {params.flags} {params.extra_flags}"
        "  {input.read1} {input.read2} 2> {log.log1}"
        " | samtools collate -u -O -"
        " | samtools fixmate -m - - 2> {log.log4}"
        " | ngs-pl filter-reads-region -o {output.bam} --outstat {log.log2} {params.mode} {params.regions} 2> {log.log3}"


# EOF
