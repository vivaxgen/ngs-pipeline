
# NOTE: this file should only define functions and variables that are used in
# the multi-sample flow (msf) rules.
# DO NOT define any rules here.

# include the global parameters
include: pkg("ngs_pipeline::global_params.smk")

# include the utilities
include: pkg("ngs_pipeline::helper/utilities.smk")

# specific for multi-sample flow, we need to know all the sample names and their indexes
# to be able to expand the input files for each sample
include: pkg("ngs_pipeline::msf/params.smk")

# specific for multi-sample flow, we need to know how to prepare the sample files (reads)
# in each sample directory.
include: pkg("ngs_pipeline::msf/prepare_sample_files.smk")

# prepare sample-related parameters


def get_sample(w):
    # return the sample name from the sample wildcard
    # check if wildcard contains sample, if not raise error
    if not hasattr(w, 'sample'):
        raise ValueError("wildcard does not contain sample")
    return w.sample


def get_indexes(w):
    # return the list of indexes for the given sample
    # check if sample is in the list of samples, if not raise error
    if w.sample not in read_files.samples():
        raise ValueError(f"sample {w.sample} not found in the list of samples")
    return read_files.get_indexes(w.sample)


def expand_index(w, pattern):
    # expand the given pattern with the list of indexes for the given sample
    # pattern should contain {idx} wildcard to be replaced by the index
    idxs = get_indexes(w)
    return expand(pattern, idx=idxs)


def add_prefix(path):
    # for multi sample flow, we need to add prefix to the path because we are running
    # outside individual sample directory. XXX: need further thinking
    return [f"{prefix}/{sample}/{path}"]


# for multi sample flow, we set sp pathvar (sample prefix) to
# {outdir}/samples/{{sample}}/

prefix = f"{outdir}/"

SP = "{pfx}/{sample}/"
pathvars:
    sp = SP


def expand_sample_index(w, template_pattern):
    # expand the given pattern with the list of indexes for the given sample
    # pattern should contain {idx} wildcard to be replaced by the index

    sample_idxs = read_files.get_indexes(w.sample)
    
    # 3. Use the Python variable here instead of the non-existent 'pathvars.sp'
    full_template = SP + template_pattern
    
    result = expand(full_template, 
                    pfx=w.pfx, 
                    sample=w.sample, 
                    idx=sample_idxs)
    return result


include: pkg("ngs_pipeline::helper/funcs.smk")


def get_mapped_bam_file():
    """ return the mapped bam file for the given sample and index """
    return replace_sp("<sp>maps/{sample}-{idx}.bam")


def get_merge_input_bam_files(w):
    """ return the list of input bam files to be merged """
    #return expand('{{pfx}}/{{sample}}/maps/mapped-final-{idx}.bam', idx=read_files.get_indexes(w.sample))
    #return expand_sample_index(w, "maps/mapped-final-{idx}.bam")
    return expand_sp('<sp>maps/mapped-final-{idx}.bam')


def get_final_file(w):
    if complete_region:
        return f'<sp>gvcf/{w.sample}-{complete_region}.g.vcf.gz'
    return [f"<sp>gvcf/{w.sample}-{reg}.g.vcf.gz" for reg in REGIONS]


# EOF
