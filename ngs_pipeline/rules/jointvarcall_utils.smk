# jointvarcall_utils.smk - ngs-pipeline rules
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2024 Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

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


# region definition

regpart = RegPartition(PARTIALS)


# final output of this workflow

def get_final_file(w):
    if complete_region:
        return f'{destdir}/vcfs/{complete_region}.vcf.gz'
    return [f"{destdir}/vcfs/{reg}.vcf.gz" for reg in REGIONS]


def get_bam_files():
    # traversing on all source directories
    return [f'{a_dir}/maps/mapped-final.bam' for (a_dir, s) in zip(SAMPLE_DIRS, SAMPLES)]


rule concat_split_vcfs:
    # this rule concatenate split-based VCF files into correspoding region
    # or chromosome-based VCF files
    threads: 2
    input:
        regpart.get_all_region_vcf
    output:
        f"{destdir}/vcfs/{{reg}}.vcf.gz"
    log:
        f"{destdir}/logs/bcftools-concat-{{reg}}.log"
    shell:
        "bcftools concat -o {output} {input} 2> {log}"


rule concat_region_vcfs:
    # this rule concatenate all region/chromosome-based VCF files
    # into single VCF file
    threads: 2
    input:
        get_final_file
    output:
        f'{destdir}/concatenated.vcf.gz'
    log:
        f'{destdir}/logs/bcftools-concat.log'
    shell:
        'bcftools concat -o {output} {input} 2> {log}'

rule csq_vcf:
    input:
        vcf = "{fn}.vcf.gz",
        tbi = "{fn}.vcf.gz.tbi",
    output:
        vcf = "{fn}-bcsq.vcf.gz",
        tbi = "{fn}-bcsq.vcf.gz.tbi"
    shell:
        "bcftools csq -l -o {output.vcf} -f {refseq} --gff {gff_file} {input.vcf}"
        " && bcftools index -t {output.vcf}"


rule concatenated_vcf:
    input:
        f'{destdir}/concatenated.vcf.gz'

# EOF
