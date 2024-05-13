# jointvarcall_split_gatk.smk - ngs-pipeline rules
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2024 Hidayat Trimarsanto <trimarsanto@gmail.com>"
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

# additional settings and parameters

# GenomicsDBImport
gdbi_flags = config.get('gdbi_flags', '')
gdbi_extra_flags = config.get('gdbi_extra_flags', '')

# GenotypeGVCFs
ggvcf_flags = config.get('ggvcf_flags', '-stand-call-conf 10 -new-qual')
ggvcf_extra_flags = config.get('ggvcf_extra_flags', '')


def get_gvcf_files(region):
    # traversing on all source directories
    return [f'{s_dir}/gvcf/{s}-{region}.g.vcf.gz' for (s_dir, s) in zip(SAMPLE_DIRS, SAMPLES)]


# list of rules

rule all:
    localrule: True
    input:
        get_final_file


# get the list of all gvcfs
rule prepare_gvcf_list:
    threads: 2
    input:
        lambda w: get_gvcf_files(w.reg)
    output:
        f"{destdir}/maps/{{reg}}.tsv"
    run:
        # write {input} to tab-delimited map file: sample\tgvcf_path
        if len(SAMPLES) != len(input):
            raise ValueError('SAMPLES and input files do not have identical length')
        with open(output[0], 'w') as f_out:
            for s, a_file in zip(SAMPLES, input):
                f_out.write(f'{s}\t{a_file}\n')


rule combine_gvcf:
    threads: thread_allocations.get('combine_gvcf', 5)
    input:
        f"{destdir}/maps/{{reg}}.tsv"
    output:
        directory(f"{destdir}/dbs/{regpart.notation}")
    log:
        f"{destdir}/logs/genomicsdbimport-{regpart.notation}.log"
    params:
        #reg = get_interval,
        reg = regpart.get_interval,
        threads = lambda w, threads: max(threads - 1, 2),
        flags = gdbi_flags,
        extra_flags = gdbi_extra_flags,
    shell:
        "gatk GenomicsDBImport --reader-threads {params.threads} "
        "--genomicsdb-workspace-path {output} --sample-name-map {input} "
        "{gdbi_flags} {gdbi_extra_flags} "
        "{params.reg} 2> {log}"


rule jointvarcall_gatk:
    threads: 3
    input:
        f"{destdir}/dbs/{regpart.notation}"
    output:
        regpart.region_vcf
    params:
        variantfile = f'--variant {variant_file}' if variant_file else '',
        reg = regpart.get_interval,
        flags = ggvcf_flags,
        extra_flags = ggvcf_extra_flags
    log:
        f"{destdir}/logs/genotypegvcfs-{regpart.notation}.log"
    shell:
        "gatk GenotypeGVCFs {params.flags} -R {refseq} -V gendb://{input} "
        "-O {output} {params.variantfile} {params.reg} {params.extra_flags} 2> {log}"


# order of the rules for {reg}.vcf.gz in case of non-split joint calling
ruleorder: jointvarcall_gatk > concat_split_vcfs


# EOF
