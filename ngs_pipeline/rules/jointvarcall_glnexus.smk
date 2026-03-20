# jointvarcall_freebayes.smk - ngs-pipeline rules
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2023-2024 Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

# jointvarcall_utils.smk will provide the following:
# variables:
#   srcdirs
#   destdir
#   SAMPLES
#   SAMPLE_DIRS
#   regpart
# functions:
#   get_final_file
#   get_bam_files
# rules:
#   concat_split_vcfs
#   concat_region_vcfs


include: "jointvarcall_utils.smk"


# list of rules

rule all:
    input:
        # get_final_file is defined in jointvarcall_utils.smk
        get_final_file


def get_gvcf_files(region):
    # traversing on all source directories
    return [f'{s_dir}/gvcf/{s}-{region}.g.vcf.gz' for (s_dir, s) in zip(SAMPLE_DIRS, SAMPLES)]

# get the list of all gvcfs
rule prepare_gvcf_list:
    localrule: True
    input:
        lambda w: get_gvcf_files(w.reg)
    output:
        f"{destdir}/maps/{{reg}}.tsv"
    run:
        # write {input} to tab-delimited map file: sample\tgvcf_path
        if len(SAMPLES) != len(input):
            raise ValueError(
                f'SAMPLES ({len(SAMPLES)}) and input files ({len(input)}) do not have identical length'
            )
        with open(output[0], 'w') as f_out:
            for s, a_file in zip(SAMPLES, input):
                f_out.write(f'{a_file}\n')


rule jointvarcall_glnexus:
    # glnexus can work with multiple threads
    threads: 16
    input:
        f"{destdir}/maps/{regpart.notation}.tsv"
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
    shell:
        "glnexus_cli --config DeepVariant"
        "  --threads {threads}"
        "  --list {input}"
        "| bcftools view -o {output}"


# order of the rules for {reg}.vcf.gz in case of non-split joint calling
ruleorder: jointvarcall_glnexus > concat_split_vcfs

# EOF
