# params and functions to initialize the single-sample variant calling workflow

# prepare necessary global parameters
from ngs_pipeline.rules import inc

include: inc("global_params.smk")

include: inc("utilities.smk")
# prepare sample-related parameters

sample = pathlib.Path.cwd().name
IDXS, = glob_wildcards('reads/raw-{idx}_R1.fastq.gz')

# for single sample flow, we set sp (sample prefix) pathvar to empty string
SP = ""
pathvars:
    sp = ""

prefix = ""

include: inc("handler_funcs.smk")


def get_sample(w):
    return sample


def get_indexes(w):
    return IDXS

def add_prefix(path):
    # for single sample flow, we don't need to add prefix to the path because we are running
    # in the sample directory, so just return the path as is
    return [path]

def get_mapped_bam_file():
    """ return the mapped bam file for the given sample and index """
    return "maps/{sample}-{{idx}}.bam".format(sample=sample)

def get_final_bam_files(w):
    """ return the final bam file for further processing """
    return ("maps/mapped-dedup-{idx}.bam" if deduplicate else "maps/mapped-filtered-{idx}.bam")

def get_merge_input_bam_files(w):
    """ return the list of input bam files to be merged """
    return expand('maps/mapped-final-{idx}.bam', idx=IDXS)

def get_final_file(w):
    if complete_region:
        return f'gvcf/{sample}-{complete_region}.g.vcf.gz'
    return [f"gvcf/{sample}-{reg}.g.vcf.gz" for reg in REGIONS]


# EOF
