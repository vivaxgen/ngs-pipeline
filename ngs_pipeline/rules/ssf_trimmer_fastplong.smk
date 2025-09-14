
fastplong_cut_tail_window_size = 10
fastplong_cut_tail_mean_quality = min_avg_quality
fastplong_trim_front = 20
fastplong_trim_tail = 20

rule reads_trimming_fastplong:
    """ use fastplong to trim reads
    """
    threads: thread_allocations.get('trimming', 8)
    input:
        read = "reads/raw-{idx}_R0.fastq.gz",
    output:
        trimmed = temp("trimmed-reads/trimmed-{idx}_R0.fastq.gz"),
    log:
        log1 = "logs/reads_trimming-{idx}.log",
        log2 = "logs/fastp-{idx}.json",
        log3 = "logs/fastp-{idx}.html"
    params:
        sample = sample,
        length_arg = f'--length_limit {maxlen}' if maxlen > 0 else '-L',
        minlen_arg = f'--length_required {minlen}' if minlen > 0 else '',
        qual_arg = f"-q {min_read_quality}" if min_read_quality > 0 else '-Q',
        cut_tail = f"--cut_tail --cut_tail_window_size {fastplong_cut_tail_window_size} --cut_tail_mean_quality {fastplong_cut_tail_mean_quality}",
        trim = f"--trim_front {fastplong_trim_front} --trim_tail {fastplong_trim_tail}",
        adapter_arg = "--disable_adapter_trimming"
    shell:
        "fastplong -w 8 {params.adapter_arg} {params.length_arg} {params.minlen_arg}"
        "  {params.qual_arg} {params.cut_tail} {params.trim}"
        "  -o {output.trimmed} -i {input.read}"
        "  -j {log.log2} -h {log.log3} > {log.log1}"


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
