
# prepare global params

include: "global_params.smk"
include: "utilities.smk"


# source directories would be provided using config=dict() args of snakemake()

srcdirs = set(config['srcdirs'])
destdir = config.get('destdir', 'vcfs')

# get all samples and sample directories

SAMPLES = []
SAMPLE_DIRS = []
for a_dir in srcdirs:
    S, = glob_wildcards(a_dir + '/{sample,[\\w-]+}')
    # filter for non-sample directories/files
    S = [s for s in S if s != 'config.yaml']
    SAMPLES += S
    SAMPLE_DIRS += [f'{a_dir}/{s}' for s in S]


# additional settings and parameters

variant_file = config.get('variant_file', None)
interval_file = config.get('interval_file', None)
ggvcf_extra_flags = config.get('ggvcf_extra_flags', '')


# final output of this workflow

def get_final_file(w):
    return [f"{destdir}/{reg}.vcf.gz.csi" for reg in REGIONS]


def get_gvcf_files(region):
    # traversing on all source directories
    return [f'{s_dir}/gvcf/{s}-{region}.g.vcf.gz' for (s_dir, s) in zip(SAMPLE_DIRS, SAMPLES)]


# define local rules

localrules: all, prepare_gvcf_files


# list of rules

rule all:
    input:
        get_final_file


# get the list of all gvcfs
rule prepare_gvcf_files:
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
    threads: 8
    input:
        f"{destdir}/maps/{{reg}}.tsv"
    output:
        directory(f"{destdir}/dbs/{{reg}}")
    shell:
        "gatk GenomicsDBImport --reader-threads 5 --genomicsdb-workspace-path {output} "
        "--sample-name-map {input} -L {wildcards.reg}"


rule jointvarcall_gatk:
    threads: 3
    input:
        f"{destdir}/dbs/{{reg}}"
    output:
        f"{destdir}/{{reg}}.vcf.gz"
    params:
        variantfile = f'--variant {variant_file}' if variant_file else '',
        intervalfile = f'--intervals {interval_file}' if interval_file else '',
    shell:
        "gatk GenotypeGVCFs -stand-call-conf 10 -new-qual -R {refseq} -V gendb://{input} "
        "-O {output} {params.variantfile} {params.intervalfile} {ggvcf_extra_flags}"

# EOF
