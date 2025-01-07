



rule reads_trimming:
    threads: 8
    input:
        read1 = "{pfx}/{sample}/trimmed-reads/dedup-{idx}_R1.fastq.gz" if optdedup else "reads/raw-{idx}_R1.fastq.gz",
        read2 = "{pfx}/{sample}/trimmed-reads/dedup-{idx}_R2.fastq.gz" if optdedup else "reads/raw-{idx}_R2.fastq.gz"
    output:
        trimmed1 = temp("{pfx}/{sample}/trimmed-reads/trimmed-{idx}_R1.fastq.gz"),
        trimmed2 = temp("{pfx}/{sample}/trimmed-reads/trimmed-{idx}_R2.fastq.gz")
    log: "logs/reads_trimming-{idx}.log"
    params:
        nextseq_arg = '--nextseq-trim 20' if is_nextseq_or_novaseq() else '',
        length_arg = f'--length {maxlen}' if maxlen > 0 else '',
        minlen_arg = f'-m {minlen}',
        qual_arg = f"-q {config['min_read_qual']}",
        adapter_arg = '' if config['libprep'] == 'generic' else adapter_arguments(config['libprep'].lower())
    shell:
        "cutadapt {params.nextseq_arg} -j {threads}"
        "  {params.length_arg} {params.minlen_arg} {params.qual_arg}"
        "  -O 3 {params.adapter_arg} -o {output.trimmed1} -p {output.trimmed2}"
        "  {input.read1} {input.read2} > {log}"

# EOF
