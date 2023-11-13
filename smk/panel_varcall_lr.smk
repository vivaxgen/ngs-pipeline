
# need to read targeted variant calling global params

import os
import pathlib

from ngsutils import cexit

include: "utilities.smk"

ngs_pipeline_basedir = os.environ['NGS_PIPELINE_BASE']
ngsenv_basedir = os.environ['NGSENV_BASEDIR']

refmap = ngsenv_basedir + '/' + config['refmap_file']
refseq = ngsenv_basedir + '/' + config['refseq_file']
target_regions = ngsenv_basedir + '/' + config['target_regions']
target_variants = ngsenv_basedir + '/' + config['target_variants']
variant_info = ngsenv_basedir + '/' + config.get('variant_info', 'NOFILE')

infiles = config['infiles']
outdir = config['outdir']
freebayes_extra_flags = config.get('freebayes_extra_flags', '')

min_read_qual = config['min_read_qual']
max_read_len = config['max_read_len']
min_read_len = config['min_read_len']
headcrop = config.get('headcrop', 0)
tailcrop = config.get('tailcrop', 0)

read_files = {}
err_files = []
for infile in infiles:
    if not infile.endswith('.fastq.gz'):
        err_files.append(
            f'Input file: {infile} is not a compressed FASTQ file (fastq.gz).'
        )
        continue
    path = pathlib.Path(infile)
    if not path.exists():
        err_files.append(
            f'Input file: {infile} does not exist. Please check your path.'
        )
        continue
    sample = path.name.removesuffix('.fastq.gz')
    read_files[sample] = infile

if any(err_files):
    cexit('ERROR: invalid input files: \n' +
          '\n'.join(f'  {errmsg}' for errmsg in err_files))


def get_read_file(w):
    return read_files[w.sample]


def get_output_file(w):
    return [f'{outdir}/{sample}/vcfs/variants.vcf.gz.csi' for sample in read_files.keys()]


wildcard_constraints:
    sample = '[\.\\w-]+'


rule all:
    localrule: True
    input:
        get_output_file


rule merged_report:
    localrule: True
    input:
        f'{outdir}/merged_report.tsv'


rule link_reads:
    localrule: True
    input:
        get_read_file
    output:
        f"{outdir}/{{sample}}/reads/raw.fastq.gz"
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
    params:
        headcrop = f'--headcrop {headcrop}' if headcrop else '',
        tailcrop = f'--tailcrop {tailcrop}' if tailcrop else '',
    shell:
        "gunzip -c {input} | chopper -t {threads} -q {min_read_qual} --minlength {min_read_len} --maxlength {max_read_len} "
        "{params.headcrop} {params.tailcrop} "
        "| gzip -c > {output}"


rule mapping:
    threads: 8
    input:
        "{pfx}/{sample}/trimmed_reads/trimmed.fastq.gz"
    output:
        "{pfx}/{sample}/maps/sorted.bam"
    params:
        rg = lambda w: f"-R @RG\\\\tID:{w.sample}\\\\tSM:{w.sample}\\\\tLB:LIB-{w.sample}",
        threads = lambda wildcards, threads: threads - 1,
    shell:
        "minimap2 -t {params.threads} -a {refmap} {params.rg} {input} "
        "|  samtools sort -@4 -o {output} "


rule freebayes:
    threads: 2
    input:
        bam = "{pfx}/maps/sorted.bam",
        idx = "{pfx}/maps/sorted.bam.bai"
    output:
        "{pfx}/vcfs/variants.vcf.gz"
    shell:
        "freebayes -f {refseq} --target {target_variants} --haplotype-length 0 --report-monomorphic "
        "--min-base-quality {min_read_qual} {freebayes_extra_flags} {input.bam} "
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


rule gen_report:
    threads: 1
    input:
        '{pfx}/vcfs/variants.vcf.gz'
    output:
        '{pfx}/genetic_report.tsv'
    shell:
        "ngs-pl generate-variant-report -o {output} --infofile {variant_info} {input}"

rule merge_report:
    localrule: True
    input:
        expand(f'{outdir}/{{sample}}/genetic_report.tsv', sample=list(read_files.keys()))
    output:
        f'{outdir}/merged_report.tsv'
    run:
        import pandas as pd

        dfs = []
        for infile in input:
            try:
                df = pd.read_table(infile)
                dfs.append(df)
            except:
                raise RuntimeError(f'ERROR parsing {infile}')

        concatenated_df = pd.concat(dfs)
        concatenated_df.to_csv(output[0], index=False, sep='\t')
        
# EOF
