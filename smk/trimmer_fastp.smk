
adapters = {
    'nextera': ['CTGTCTCTTATACACATCT', 'CTGTCTCTTATACACATCT'],
    'truseq': ['AGATCGGAAGAGCACACGTCTGAACTCCAGTCA', 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'],
    'truseqstranded': ['AGATCGGAAGAGCACACGTCTGAACTCCAGTCA', 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'],
    'ultraii': ['AGATCGGAAGAGCACACGTCTGAACTCCAGTCA', 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'],
    'ultraiifs': ['AGATCGGAAGAGCACACGTCTGAACTCCAGTCA', 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'],
    'dnaprep': ['CTGTCTCTTATACACATCT+ATGTGTATAAGAGACA', 'CTGTCTCTTATACACATCT+ATGTGTATAAGAGACA'],
}
minlen = int(config['minlen']) if 'minlen' in config else int(config['read_length'] / 3)
maxlen = int(config['maxlen']) if 'maxlen' in config else 0


def adapter_arguments(libprep):
    arguments = []
    adapter_1, adapter_2 = adapters[libprep]
    arguments += ['-a %s' % x for x in adapter_1.split('+')]
    arguments += ['-A %s' % x for x in adapter_2.split('+')]
    return ' '.join(arguments)


def is_nextseq_or_novaseq():
    return config['instrument'].lower().startswith('nextseq') or config['instrument'].lower().startswith('novaseq')


optdedup = config.get('optical_dedup', False)


rule reads_trimming:
    threads: 8
    input:
        read1 = "trimmed-reads/dedup-{idx}_R1.fastq.gz" if optdedup else "reads/raw-{idx}_R1.fastq.gz",
        read2 = "trimmed-reads/dedup-{idx}_R2.fastq.gz" if optdedup else "reads/raw-{idx}_R2.fastq.gz"
    output:
        trimmed1 = temp("trimmed-reads/trimmed-{idx}_R1.fastq.gz"),
        trimmed2 = temp("trimmed-reads/trimmed-{idx}_R2.fastq.gz")
    log:
        log1 = "logs/reads_trimming-{idx}.log",
        log2 = "logs/fastp-{idx}.json",
        log3 = "logs/fastp-{idx}.html"
    params:
        nextseq_arg = '--trim_poly_g' if is_nextseq_or_novaseq() else '',
        length_arg = f'--length_limit {maxlen}' if maxlen > 0 else '',
        minlen_arg = f'--length_required {minlen}',
        qual_arg = f"-q {config['min_read_qual']}",
        adapter_arg = '' if config['libprep'] == 'generic' else adapter_arguments(config['libprep'].lower()),
        correction_arg = '-c' if config.get('correction', False) else ''
    shell:
        "fastp -w 8 {params.correction_arg} {params.nextseq_arg} {params.length_arg} {params.minlen_arg} {params.qual_arg} -o {output.trimmed1} -O {output.trimmed2} -i {input.read1} -I {input.read2} -j {log.log2} -h {log.log3} > {log.log1}"

rule trimming_stat:
    localrule: True
    input:
        "logs/fastp-{idx}.json"
    output:
        "logs/trimming_stat-{idx}.json"
    run:
        import json

        fastp_d = json.load(open(input[0]))
        d = dict(
            original_reads=fastp_d['summary']['before_filtering']['total_reads'],
            filtered_reads=fastp_d['summary']['after_filtering']['total_reads']
        )

        json.dump(d, open(output[0], 'w'))

# EOF
