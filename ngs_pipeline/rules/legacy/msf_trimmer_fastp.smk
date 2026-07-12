

def is_nextseq_or_novaseq():
    return config['instrument'].lower().startswith('nextseq') or config['instrument'].lower().startswith('novaseq')


optdedup = config.get('optical_dedup', False)


rule reads_trimming:
    threads: thread_allocations.get('trimming', 8)
    input:
        read1 = "{pfx}/trimmed-reads/dedup-{idx}_R1.fastq.gz" if optdedup else "{pfx}/reads/raw-{idx}_R1.fastq.gz",
        read2 = "{pfx}/trimmed-reads/dedup-{idx}_R2.fastq.gz" if optdedup else "{pfx}/reads/raw-{idx}_R2.fastq.gz"
    output:
        trimmed1 = temp("{pfx}/trimmed-reads/trimmed-{idx}_R1.fastq.gz"),
        trimmed2 = temp("{pfx}/trimmed-reads/trimmed-{idx}_R2.fastq.gz")
    log:
        log1 = "{pfx}/logs/reads_trimming-{idx}.log",
        log2 = "{pfx}/logs/fastp-{idx}.json",
        log3 = "{pfx}/logs/fastp-{idx}.html"
    params:
        nextseq_arg = '--trim_poly_g' if is_nextseq_or_novaseq() else '',
        length_arg = f'--length_limit {maxlen}' if maxlen > 0 else '-L',
        minlen_arg = f'--length_required {minlen}' if minlen > 0 else '',
        qual_arg = f"-q {min_read_quality}" if min_read_quality > 0 else '-Q',
        avg_qual_arg = f"-e {min_avg_quality}" if min_avg_quality > 0 else '',
        adapter_arg = '' if config['libprep'] == 'generic' else '-A',
        correction_arg = '-c' if config.get('correction', False) else ''
    shell:
        "fastp -w 8 {params.correction_arg} {params.nextseq_arg} {params.adapter_arg} {params.length_arg} {params.minlen_arg} {params.qual_arg} -o {output.trimmed1} -O {output.trimmed2} -i {input.read1} -I {input.read2} -j {log.log2} -h {log.log3} > {log.log1}"


rule trimming_stat:
    localrule: True
    input:
        "{pfx}/logs/fastp-{idx}.json"
    output:
        "{pfx}/logs/trimming_stat-{idx}.json"
    run:
        import json

        fastp_d = json.load(open(input[0]))
        d = dict(
            original_reads=fastp_d['summary']['before_filtering']['total_reads'],
            filtered_reads=fastp_d['summary']['after_filtering']['total_reads']
        )

        json.dump(d, open(output[0], 'w'))

# EOF
