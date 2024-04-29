# jointvarcall_freebayes.smk - ngs-pipeline rules
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2023-2024 Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

# jointvarcall_utils.smk will provide the following:
# variables:
#   srcdirs
#   destdir
#   SAMPLE_DIRS
#   regpart
# functions:
#   get_final_file
#   get_bam_files
# rules:
#   concat_split_vcfs
#   concat_region_vcfs


include: "jointvarcall_utils.smk"

# flags is for mandatory arguments
freebayes_flags = config.get('freebayes_flags', '--min-alternate-count 2')

# extra_flags is user-modifiable arguments
freebayes_extra_flags = config.get('freebayes_extra_flags', '')

# freebayes argument is --region chrom:start-end or --targets bed_file
# so we need to setup the regpart using
# set_arg_name(region_argument_name, bedfile_argument_name)
regpart.set_arg_name('--region', '--targets')


# list of rules

rule all:
    input:
        # get_final_file is defined in jointvarcall_utils.smk
        get_final_file


# get the list of all bams
rule prepare_bam_list:
    localrule: True
    input:
        # get_bam_files() is defined in jointvarcall_utils.smk
        lambda w: get_bam_files()
    output:
        f"{destdir}/meta/bam_list.txt"
    run:
        # check format for freebayes bam file list
        # write {input} to ??
        if len(SAMPLES) != len(input):
            raise ValueError('SAMPLES and input files do not have identical length')
        with open(output[0], 'w') as f_out:
            for s, a_file in zip(SAMPLES, input):
                f_out.write(f'{a_file}\n')


rule jointvarcall_freebayes:
    # freebayes can only work with single thread
    threads: 1
    input:
        f"{destdir}/meta/bam_list.txt"
    output:
        # regpart.region_vcf() will return either:
        #    temp(f"{destdir}/split/{{reg}}~{{idx}}.vcf.gz")
        # or
        #    f"{destdir}/vcfs/{{reg}}.vcf.gz"
        regpart.region_vcf
    params:
        # regpart.get_interval() will return either one of the following:
        # --region CHROM, --region CHROM:START-END, --targets BEDFILE_PATH
        reg = regpart.get_interval,
        flags = freebayes_flags + ' ' + freebayes_extra_flags,
    shell:
        "freebayes --fasta-reference {refseq} -L {input} --ploidy {ploidy} "
        "{params.reg} {params.flags} "
        "| bcftools view -o {output}"

# EOF
