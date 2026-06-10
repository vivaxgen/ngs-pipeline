
fastplong_cut_tail_window_size = config.get('fastplong_cut_tail_window_size', 10)
fastplong_cut_tail_mean_quality = config.get('fastplong_cut_tail_mean_quality', 1)
fastplong_trim_front = config.get('fastplong_trim_front', 20)
fastplong_trim_tail = config.get('fastplong_trim_tail', 20)


rule msf_trim_reads:
    threads: 4
    input:
        read = "{pfx}/{sample}/reads/raw-{idx}.fastq.gz"
    output:
        trimmed = "{pfx}/{sample}/trimmed-reads/trimmed-{idx}.fastq.gz"
    log:
        log1 = "{pfx}/{sample}/logs/reads_trimming-{idx}.log",
        log2 = "{pfx}/{sample}/logs/fastp-{idx}.json",
        log3 = "{pfx}/{sample}/logs/fastp-{idx}.html"
    params:
        length_arg = f'--length_limit {maxlen}' if maxlen > 0 else '-L',
        minlen_arg = f'--length_required {minlen}' if minlen > 0 else '',
        qual_arg = f"-q {min_read_qual}" if min_read_qual > 0 else '-Q',
        cut_tail = f"--cut_tail --cut_tail_window_size {fastplong_cut_tail_window_size} --cut_tail_mean_quality {fastplong_cut_tail_mean_quality}",
        trim = f"--trim_front {fastplong_trim_front} --trim_tail {fastplong_trim_tail}",
        adapter_arg = "--disable_adapter_trimming"
    shell:
        "fastplong -w {threads} {params.adapter_arg} {params.length_arg} {params.minlen_arg}"
        "  {params.qual_arg} {params.cut_tail} {params.trim}"
        "  -o {output.trimmed} -i {input.read}"
        "  -j {log.log2} -h {log.log3} > {log.log1}"

# EOF