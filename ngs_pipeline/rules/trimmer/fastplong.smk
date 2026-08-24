# SPDX-FileCopyrightText: 2024-2006 Hidayat (Anto)Trimarsanto <trimarsanto@gmail.com>
# SPDX-License-Identifier: MIT

__copyright__ = "(c) 2024-2006 Hidayat (Anto) Trimarsanto <trimarsanto@gmail.com>"
__author__ = "trimarsanto@gmail.com"
__license__ = "MIT"

fastplong_cut_tail_window_size = config.get('fastplong_cut_tail_window_size', 10)
fastplong_cut_tail_mean_quality = config.get('fastplong_cut_tail_mean_quality', 1)
fastplong_trim_front = config.get('fastplong_trim_front', 20)
fastplong_trim_tail = config.get('fastplong_trim_tail', 20)


rule reads_trimming_lr:
    threads: 4
    input:
        read = "{anypath}reads/raw-{idx}_R0.fastq.gz",
        model = "{anypath}reads/model-{idx}.txt" if ngs_platform.upper() in ['ONT'] else []
    output:
        trimmed = temp_unless("{anypath}trimmed-reads/trimmed-{idx}.fastq.gz", keep_trimmed_fastq)
    log:
        log1 = "{anypath}logs/reads_trimming-{idx}.log",
        log2 = "{anypath}logs/fastplong-{idx}.json",
        log3 = "{anypath}logs/fastplong-{idx}.html"
    params:
        length_arg = f'--length_limit {maxlen}' if maxlen > 0 else '-L',
        minlen_arg = f'--length_required {minlen}' if minlen > 0 else '',
        qual_arg = f"-q {min_read_qual}" if min_read_qual > 0 else '-Q',
        cut_tail = (f"--cut_tail --cut_tail_window_size {fastplong_cut_tail_window_size} --cut_tail_mean_quality {fastplong_cut_tail_mean_quality}"
                    if fastplong_cut_tail_window_size > 0 and fastplong_cut_tail_mean_quality > 0 else ''),
        trim = (f"--trim_front {fastplong_trim_front} --trim_tail {fastplong_trim_tail}"
                if fastplong_trim_front > 0 and fastplong_trim_tail > 0 else ''),
        adapter_arg = "--disable_adapter_trimming"
    shell:
        "fastplong -w {threads} {params.adapter_arg} {params.length_arg} {params.minlen_arg}"
        "  {params.qual_arg} {params.cut_tail} {params.trim}"
        "  -o {output.trimmed} -i {input.read}"
        "  -j {log.log2} -h {log.log3} > {log.log1}"


rule trimming_stat:
    localrule: True
    input:
        "{anypath}logs/fastplong-{idx}.json"
    output:
        "{anypath}logs/trimming_stat-{idx}.json"
    run:
        import json

        fastp_d = json.load(open(input[0]))
        d = dict(
            original_reads=fastp_d['summary']['before_filtering']['total_reads'],
            filtered_reads=fastp_d['summary']['after_filtering']['total_reads']
        )

        json.dump(d, open(output[0], 'w'))

# EOF