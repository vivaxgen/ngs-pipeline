
# need to read targeted variant calling global params

import os
import pathlib

ngs_pipeline_basedir = os.environ['NGS_PIPELINE_BASE']
ngsenv_basedir = os.environ['NGSENV_BASEDIR']

refmap = ngsenv_basedir + '/' + config['refmap_file']
refseq = ngsenv_basedir + '/' + config['refseq_file']
target_regions = ngsenv_basedir + '/' + config['target_regions']
target_variants = ngsenv_basedir + '/' + config['target_variants']

infiles = config['infiles']
outdir = config['outdir']

max_read_len = config['max_read_len']
min_read_len = config['min_read_len']

read_files = {}
for infile in infiles:
    path = pathlib.Path(infile)
    sample = path.name.removesuffix('.fastq.gz')
    read_files[sample] = infile


def get_read_file(w):
    return read_files[w.sample_id]


def get_output_file(w):
    return [f'{outdir}/{sample_id}/vcfs/variants.vcf.gz.csi' for sample_id in read_files.keys()]


wildcard_constraints:
    sample_id = '[\.\\w-]+'


rule all:
    localrule: True
    input:
        get_output_file


rule link_reads:
    localrule: True
    input:
        get_read_file
    output:
        f"{outdir}/{{sample_id}}/reads/raw.fastq.gz"
    run:
        import pathlib

        dest_file = pathlib.Path(output[0])
        src_file = pathlib.Path(input[0]).resolve()
        dest_file.symlink_to(src_file)


rule plot_qc:
    threads: 2
    input:
        "{pfx}/reads/raw.fastq.gz"
    output:
        directory("{pfx}/QC/")
    shell:
        "NanoPlot -t 2 --fastq {input} --outdir {output}"


rule trim_reads:
    threads: 4
    input:
        "{pfx}/reads/raw.fastq.gz"
    output:
        "{pfx}/trimmed_reads/trimmed.fastq.gz"
    shell:
        "gunzip -c {input} | chopper -t {threads} -q 10 --minlength {min_read_len} --maxlength {max_read_len} | gzip -c > {output}"


rule mapping:
    threads: 8
    input:
        "{pfx}/{sample_id}/trimmed_reads/trimmed.fastq.gz"
    output:
        "{pfx}/{sample_id}/maps/sorted.bam"
    params:
        rg = lambda w: f"-R @RG\\\\tID:{w.sample_id}\\\\tSM:{w.sample_id}\\\\tLB:LIB-{w.sample_id}",
    shell:
        "minimap2 -a {refmap} {params.rg} {input} "
        "|  samtools sort -@4 -o {output} "


rule freebayes:
    threads: 1
    input:
        bam = "{pfx}/maps/sorted.bam",
        idx = "{pfx}/maps/sorted.bam.bai"
    output:
        "{pfx}/vcfs/variants.vcf.gz"
    shell:
        "freebayes -f {refseq} --target {target_variants} --report-monomorphic {input.bam} "
        "| bcftools sort -o {output}"


rule filter_vcf:
    threads: 1
    input:
        "{pfx}/vcfs/raw.bayes.vcf.gz"
    output:
        "{pfx}/vcfs/hq.bayes.vcf.gz"
    shell:
        "bcftools view -o {output} -e 'QUAL<15' {input} "
        "&& sleep 1 && tabix {output}"


rule freebayes_2:
    threads: 1
    input:
        bam = "{pfx}/maps/sorted.bam",
        base = "{pfx}/vcfs/hq.bayes.vcf.gz"
    output:
        "{pfx}/vcfs/haplotypes.vcf.gz"
    shell:
        "freebayes -f {refseq} --haplotype-length 500 --target {target_regions} --haplotype-basis-alleles {input.base} {input.bam} "
        "| bgzip -c > {output}"


# utilities

rule index_csi:
    threads: 1
    input:
        "{pfx}.vcf.gz"
    output:
        "{pfx}.vcf.gz.csi"
    shell:
        "bcftools index --csi {input}"


rule index_tbi:
    threads: 1
    input:
        "{pfx}.vcf.gz"
    output:
        "{pfx}.vcf.gz.tbi"
    shell:
        "bcftools index --tbi {input}"

rule index_bai:
    threads: 1
    input:
        "{pfx}.bam"
    output:
        "{pfx}.bam.bai"
    shell:
        "samtools index {input}"


# EOF
