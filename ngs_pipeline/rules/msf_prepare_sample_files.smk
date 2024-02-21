# msf_prepare_sample_files.smk - ngs-pipeline rules
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2023, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

import pathlib

rule msf_link_reads_se:
    localrule: True
    input:
        lambda w: read_files.get_read_file(w)
    output:
        f"{outdir}/{{sample}}/reads/raw-{{idx}}.fastq.gz"
    run:

        dest_file = pathlib.Path(output[0])
        src_file = pathlib.Path(input[0]).resolve()
        dest_file.symlink_to(src_file)


rule msf_link_reads_pe:
    localrule: True
    input:
        lambda w: read_files.get_read_file(w)
    output:
        f"{outdir}/{{sample}}/reads/raw-{{idx}}_R1.fastq.gz",
        f"{outdir}/{{sample}}/reads/raw-{{idx}}_R2.fastq.gz"
    run:

        if (count := len(input)) != 2:
            raise ValueError(
                f'Trying to link in paired-end mode for sample {wildcards.sample}, '
                f'but number of input file is {count}')

        for src, dst in zip(input, output):
            dest_file = pathlib.Path(dst)
            src_file = pathlib.Path(src).resolve()
            dest_file.symlink_to(src_file)  


# EOF
