
# prepare global params

include: "global_params.smk"


# source directories would be provided using config=dict() args of snakemake()

srcdirs = config['srcdirs']
destdir = config.get('destdir', 'vcfs')

# get all samples and sample directories

SAMPLES = []
SAMPLE_DIRS = []
for a_dir in srcdirs:
    S, = glob_wildcards(a_dir + '/{sample,[\\w-]+}')
    SAMPLES += [s for s in S if s != 'config.yaml']
    SAMPLE_DIRS += [f'{a_dir}/{s}' for s in S]


# final output of this workflow

def get_final_file(w):
    return [f"{destdir}/joint-{reg}.vcf.gz" for reg in REGIONS]


def get_bam_files():
    # traversing on all source directories
    return [f'{a_dir}/maps/mapped-final.bam' for (a_dir, s) in zip(SAMPLE_DIRS, SAMPLES)]


# define local rules

localrules: all, prepare_bam_list


# list of rules

rule all:
    input:
        get_final_file


# get the list of all bams
rule prepare_bam_list:
    threads: 2
    input:
        lambda w: get_bam_files()
    output:
        f"{destdir}/bam_list.txt"
    run:
        # check format for freebayes bam file list
        # write {input} to ??
        if len(SAMPLES) != len(input):
            raise ValueError('SAMPLES and input files do not have identical length')
        with open(output[0], 'w') as f_out:
            for s, a_file in zip(SAMPLES, input):
                f_out.write(f'{a_file}\n')

rule jointvarcall_freebayes:
    threads: 3
    input:
        f"{destdir}/bam_list.txt"
    output:
        f"{destdir}/joint-{{reg}}.vcf.gz"
    shell:
        "freebayes --fasta-reference {refseq} --ploidy {ploidy} --min-alternate-count 2 --region {wildcards.reg} > {output}"
        " && sleep 1 && bcftools index {output}"

# EOF
