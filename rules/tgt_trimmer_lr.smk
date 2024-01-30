
# params:
# - min_read_qual
# - min_read_len
# - max_read_len
# - headcrop
# - tailcrop


rule plot_qc:
    threads: 2
    input:
        raw = "{pfx}/reads/raw-{idx}.fastq.gz"
    output:
        dir = directory("{pfx}/QC-{idx}/")
    log:
        err = "{pfx}/logs/nanoplot-{idx}.log"
    shell:
        "NanoPlot -t 2 --fastq {input.raw} --outdir {output} 2> {log.err}"


rule trim_reads:
    threads: 4
    input:
        raw = "{pfx}/reads/raw-{idx}.fastq.gz"
    output:
        fastq = "{pfx}/trimmed_reads/trimmed-{idx}.fastq.gz"
    params:
        headcrop = f'--headcrop {headcrop}' if headcrop else '',
        tailcrop = f'--tailcrop {tailcrop}' if tailcrop else '',
    log:
        err = "{pfx}/logs/chopper-{idx}.log"
    shell:
        "gunzip -c {input} | chopper -t {threads} -q {min_read_qual} --minlength {min_read_len} --maxlength {max_read_len} "
        "{params.headcrop} {params.tailcrop} 2> {log.err}"
        "| gzip -c > {output}"

# EOF
