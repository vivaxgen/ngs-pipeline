
# prepare global params

include: "global_params.smk"


# source directories would be provided using config=dict() args of snakemake()

srcdirs = config['srcdirs']

# get all samples and sample directories

SAMPLES = []
SAMPLE_DIRS = []
for a_dir in srcdirs:
    S, = glob_wildcards(a_dir + '/{sample,[\\w-]+}')
    SAMPLES += S
    SAMPLE_DIRS += [f'{a_dir}/{s}' for s in S]


# final output of this workflow

def get_final_file(w):
    return [f"vcfs/joint-{reg}.vcf.gz" for reg in REGIONS]


def get_gvcf_files(region):
    # traversing on all source directories
    return [f'{a_dir}/gvcf/{s}-{region}.g.vcf' for (a_dir, s) in zip(SAMPLE_DIRS, SAMPLES)]


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
        "vcfs/maps/{reg}.tsv"
    run:
        # write {input} to tab-delimited map file: sample\tgvcf_path
        if len(SAMPLES) != len(input):
            raise ValueError('SAMPLES and input files do not have identical length')
        with open(output[0], 'w') as f_out:
            for s, a_file in zip(SAMPLES, input):
                f_out.write(f'{s}\t{a_file}\n')


rule combine_gvcf:
    threads: 6
    input:
        "vcfs/maps/{reg}.tsv"
    output:
        directory("vcfs/dbs/{reg}")
    shell:
        #"touch {output}"
        "gatk GenomicsDBImport --reader-threads 5 --genomicsdb-workspace-path {output} --sample-name-map {input} -L {wildcards.reg}"


rule jointvarcall_gatk:
    threads: 3
    input:
        directory("vcfs/dbs/{reg}")
    output:
        "vcfs/joint-{reg}.vcf.gz"
    shell:
        #"touch {output}"
        "gatk GenotypeGVCFs -stand-call-conf 10 -new-qual -R {refseq} -V gendb://{input} -O {output}"

# EOF

