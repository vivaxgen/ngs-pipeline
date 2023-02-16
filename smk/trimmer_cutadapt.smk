
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


def is_nextseq():
    return config['instrument'].lower().startswith('nextseq')


optdedup = config.get('optical_dedup', False)


rule optical_dedup:
    input:
        read1 = "reads/raw-{idx}_R1.fastq.gz",
        read2 = "reads/raw-{idx}_R2.fastq.gz"
    output:
        dedup1 = "trimmed-reads/dedup-{idx}_R1.fastq.gz",
        dedup2 = "trimmed-reads/dedup-{idx}_R2.fastq.gz"
    log: "logs/optical_dedup-{idx}.log"
    shell:
        "clumpify.sh in={input.read1} in2={input.read2} out1={output.dedup1} out2={output.dedup2} dedupe optical %s 2> {log}"
        % ('spany adjacent' if is_nextseq() else '')


rule reads_trimming:
    threads: 8
    input:
        read1 = "trimmed-reads/dedup-{idx}_R1.fastq.gz" if optdedup else "reads/raw-{idx}_R1.fastq.gz",
        read2 = "trimmed-reads/dedup-{idx}_R2.fastq.gz" if optdedup else "reads/raw-{idx}_R2.fastq.gz"
    output:
        trimmed1 = temp("trimmed-reads/trimmed-{idx}_R1.fastq.gz"),
        trimmed2 = temp("trimmed-reads/trimmed-{idx}_R2.fastq.gz")
    log: "logs/reads_trimming-{idx}.log"
    params:
        nextseq_arg = '--nextseq-trim 20' if is_nextseq() else '',
        length_arg = f'--length {maxlen}' if maxlen > 0 else '',
        minlen_arg = f'-m {minlen}',
        qual_arg = f"-q {config['min_read_qual']}",
        adapter_arg = '' if config['libprep'] == 'generic' else adapter_arguments(config['libprep'].lower())
    shell:
        "cutadapt {params.nextseq_arg} -j {threads} {params.length_arg} {params.minlen_arg} {params.qual_arg} -O 3 {params.adapter_arg} -o {output.trimmed1} -p {output.trimmed2} {input.read1} {input.read2} > {log}"

# EOF
