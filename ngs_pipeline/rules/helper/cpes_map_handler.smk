# SPDX-FileCopyrightText: 2026 Hidayat Trimarsanto <trimarsanto@gmail.com>
# SPDX-License-Identifier: MIT

# this file contains rules for processing mapped BAM files based on the consecutive
# file
# prepend with @ for temporary file
# {sample}-filtered{ff,rr,rf,trans,sec,sup,unmapped}@-optdup{rm}@-clipped@-target{g6pd}@-{idx}.bam

# A Chained Path-Encoded State (CPES) is a workflow design pattern in which
# the requested output path encodes the desired sequence of transformations,
# allowing the workflow engine to infer the corresponding execution graph

# The CPES grammar for this module is:
# - clip
# - dup
# - optdup
# - target
# - filter
# 
# ie: -<stage>(options)

# Note: since the filename might contains parentheses, all arguments in shell requiring the filename
# should be quoted, e.g. '{input.bam}' or "{input.bam}" to avoid shell interpretation of special characters
# Note: the following characters are not quoted by ls: # % : @ {} ~
# alternative: {sample}-filter:ff:rf:trans:sec:sup:unmapped:maxis%300_-dup:rm_-clip:soft_-target:g6pd-{idx}.bam

from time import sleep
from ngs_pipeline.rules import inc


wildcard_constraints:
    #options = r".*",
    options = r"(?::(?:[_-]?[^:_-])*)*",
    stages = r"[^/]+",


include: inc("ngs_pipeline::helper/stats.smk")


def get_final_bam_file(w):
    """ return the final bam file for further processing """
    #return "maps/mapped-filter-dup-{idx}.bam"
    #return "maps/mapped-filter(rr+ff+rf)-{idx}.bam"
    return config.get("final_bam_file","<sp>maps/mapped-filter-dup:rm-{idx}.bam")


rule map_link:
    localrule: True
    input:
        bam = get_mapped_bam_file()
    output:
        bam = temp("<sp>maps/mapped-{idx}.bam")
    shell:
        "ln -f {input.bam} {output.bam}"


rule map_final_link:
    localrule: True
    input:
        bam = get_final_bam_file
    output:
        bam = temp("<sp>maps/mapped-final-{idx}.bam")
    shell:
        "ln -f '{input.bam}' '{output.bam}'"


def get_filter_options(w):
    """Get filter options for ngs-pl filter-reads-orientation."""
    args = []

    options = w.options.removeprefix(":")
    if options == "":
        return '--remove_FF --remove_RR --remove_RF --remove_trans --remove_secondary --remove_supplementary --remove_unmapped'

    for opt in options.split(':'):
        match opt:
            case 'ff':
                args.append('--remove_FF')
            case 'rr':
                args.append('--remove_RR')
            case 'rf':
                args.append('--remove_RF')
            case 'trans':
                args.append('--remove_trans')
            case 'sec':
                args.append('--remove_secondary')
            case 'sup':
                args.append('--remove_supplementary')
            case 'unmapped':
                args.append('--remove_unmapped')
            case _ if opt.startswith('maxis='):
                size = opt.removeprefix('maxis=').removesuffix('bp')
                args.extend(['--max-insert-size', size])
            case _:
                pass

    return ' '.join(args)


rule bam_filtering:
    # note: input file should be name-sorted or name-collated (ie. directly from mapper)
    # output will be coordinate-sorted and indexed
    input:
        bam = "<sp>maps/{stages}-{idx}.bam",
        # no need for index file because the operation is not coordinate-based
    output:
        bam = "<sp>maps/{stages}-filter{options}-{idx}.bam",
        idx = "<sp>maps/{stages}-filter{options}-{idx}.bam.bai",
    log:
        log1 = "<sp>logs/{stages}-filter{options}-{idx}.log",
        log2 = "<sp>logs/{stages}-filter{options}-{idx}.json",
    params:
        sample = get_sample,
        args = get_filter_options,
    shell:
        "ngs-pl filter-reads-orientation --outstat '{log.log2}' {params.args} '{input.bam}' 2> '{log.log1}' "
        " | samtools sort -@4 -o '{output.bam}'"
        " && sleep 1"
        " && samtools index '{output.bam}'"


use rule bam_filtering as bam_filtering_tmp with:
    output:
        bam = temp("<sp>maps/{stages}-filter{options}_-{idx}.bam"),
        idx = temp("<sp>maps/{stages}-filter{options}_-{idx}.bam.bai")


def get_dup_options(w):
    """ get duplication options to samtools markdup """

    options = w.options.removeprefix(":") # remove leading colon
    if options == "":
        return "-r"

    args = []
    if 'rm' in options:
        args.append('-r')
    return ' '.join(args)


rule bam_deduplicating:
    input:
        bam = "<sp>maps/{stages}-{idx}.bam",
        idx = "<sp>maps/{stages}-{idx}.bam.bai",
    output:
        bam = "<sp>maps/{stages}-dup{options}-{idx}.bam",
        idx = "<sp>maps/{stages}-dup{options}-{idx}.bam.bai",
    log:
        log1 = "<sp>logs/{stages}-dup{options}-{idx}.log",
        log2 = "<sp>logs/{stages}-dup{options}-{idx}.json",
    params:
        sample = get_sample,
        args = get_dup_options,
    shell:
        "samtools markdup {params.args} --json -f '{log.log2}' '{input.bam}' '{output.bam}' 2> '{log.log1}'"
        " && sleep 1"
        " && samtools index '{output.bam}'"


use rule bam_deduplicating as bam_deduplicating_tmp with:
    output:
        bam = temp("<sp>maps/{stages}-dup{options}_-{idx}.bam"),
        idx = temp("<sp>maps/{stages}-dup{options}_-{idx}.bam.bai")


def get_target_options(w):
    """ get target options to samtools view """
    options = w.options.removeprefix(":")
    if options == "":
        return f"-L {targetregion_file}"
    args = []
    return ' '.join(args)


rule bam_targeting:
    input:
        bam = "<sp>maps/{stages}-{idx}.bam",
        idx = "<sp>maps/{stages}-{idx}.bam.bai",
    output:
        bam = "<sp>maps/{stages}-target{options}-{idx}.bam",
        idx = "<sp>maps/{stages}-target{options}-{idx}.bam.bai",
    log:
        log1 = "<sp>logs/{stages}-target{options}-{idx}.log",
        log2 = "<sp>logs/{stages}-target{options}-{idx}.json",
    params:
        sample = get_sample,
        args = get_target_options,
    shell:
        "samtools view -@ {threads} -o {output.bam} {params.args} {input.bam}"
        " && sleep 1"
        " && samtools index '{output.bam}'"


use rule bam_targeting as bam_targeting_tmp with:
    output:
        bam = temp("<sp>maps/{stages}-target{options}_-{idx}.bam"),
        idx = temp("<sp>maps/{stages}-target{options}_-{idx}.bam.bai")
    

def get_clip_options(w):
    """ get clipping options to samtools view """
    options = w.options.removeprefix(":")
    args = []
    if 'soft' in options:
        args.append('-h')
    return ' '.join(args)


rule bam_clipping:
    input:
        bam = "<sp>maps/{stages}-{idx}.bam",
        idx = "<sp>maps/{stages}-{idx}.bam.bai",
    output:
        bam = "<sp>maps/{stages}-clip{options}-{idx}.bam",
        idx = "<sp>maps/{stages}-clip{options}-{idx}.bam.bai",
    params:
        args = get_clip_options,
    shell:
        "samtools view {params.args} '{input.bam}' | samtools sort -@4 -o '{output.bam}' && samtools index '{output.bam}'"


rule bam_sorting:
    # note: options is pos or name
    input:
        bam = "<sp>maps/{stages}-{idx}.bam",
        idx = "<sp>maps/{stages}-{idx}.bam.bai",
    output:
        bam = "<sp>maps/{stages}-sort{options}-{idx}.bam",
        idx = "<sp>maps/{stages}-sort{options}-{idx}.bam.bai",
    shell:
        "samtools sort -@4 -o '{output.bam}' '{input}' && samtools index '{output.bam}'"


rule bam_collate:
    input:
        bam = "<sp>maps/{stages}-{idx}.bam",
        idx = "<sp>maps/{stages}-{idx}.bam.bai",
    output:
        bam = "<sp>maps/{stages}-collate{options}-{idx}.bam",
        idx = "<sp>maps/{stages}-collate{options}-{idx}.bam.bai",
    shell:
        "samtools collate -o '{output.bam}' '{input.bam}'"


rule map_merge_final:
    # this rule merges dedup input BAM
    threads: thread_allocations.get('map_merging', 4)
    input:
        #expand('<sp>maps/mapped-final-{idx}.bam', idx=IDXS)
        #get_merge_input_bam_files
        expand_sp('<sp>maps/mapped-final-{idx}.bam')
    output:
        temp_unless('<sp>maps/mapped-final.bam', keep_final_bam)
    params:
        sample = get_sample,
    run:
        if len(input) > 1:
            shell('samtools merge -@4 {output} {input}')
        else:
            # use hard link since input will be removed
            shell('ln {input} {output}')
        # to make index file newer by 1-sec to bam file
        sleep(2)
        shell('samtools index {output}')


rule collect_mapping_stats:
    # XXX: this one need fixing
    # probably use such:
    # ngs-pl collect-stats -o {output} --mindepth {min_depth} 
    threads: 1
    input:
        trims = expand_sp('<sp>logs/trimming_stat-{idx}.json'),
        maps = expand_sp('<sp>logs/{sample}-{idx}.stats.txt'),
        depths = expand_sp('<sp>logs/mapped-final-{idx}.depths.txt.gz'),
    params:
        sample = get_sample,
        trimmed = lambda wildcards, input: '--trimmed ' + ' --trimmed '.join(input.trims),
        mapped = lambda wildcards, input: '--mapped ' + ' --mapped '.join(input.maps),
        depthed = lambda wildcards, input: '--depth ' + ' --depth '.join(input.depths),
    output:
        '<sp>logs/stats.tsv_xxx'
    shell:
        "touch {output}"

ruleorder: bam_filtering_tmp > bam_filtering > index_bai


# the following is very experimental and might not work at all

import re

def parse_pipeline_logs(w):
    """
    Recursively parses the final BAM filename backwards 
    to extract and reconstruct all intermediate log file paths.
    """
    # 1. Resolve the actual final filename string from config/wildcards
    # Using getattr handles cases where 'w' is a wildcards object or standard dict
    #sample = get_sample(w)
    #idx_val = getattr(w, 'idx', '0')
    final_files = expand_sp(get_final_bam_file(w))(w)

    print(final_files)

    collected_logs = {}
    file_pattern = re.compile(r"/([^/]+?)-(\d+)\.bam$")

    # Regex to capture the trailing stage pattern: -stage{options} or -stage{options}_
    # Matches strings like '-filter:ff', '-dup:rm_-', etc.
    stage_pattern = re.compile(r"-([a-zA-Z0-9]+)((?::(?:[_-]?[^:_-])*)*)(_?)$")

    for final_file in final_files:
    
        # Strip paths to focus only on the naming chain and get the current idx
        file_match = file_pattern.search(final_file)
        if file_match:
            base_name = file_match.group(1)
            idx = file_match.group(2)
        else:
            raise ValueError(f"Unexpected final BAM filename format: {final_file}")
                
        current_string = base_name

        # 2. Walk backwards up the compositional chain
        while True:
            match = stage_pattern.search(current_string)
            if not match:
                break  # Reached the root file (e.g., 'mapped')
                
            stage = match.group(1)     # e.g., 'filter', 'dup', 'target'
            options = match.group(2)   # e.g., ':ff', ':rm'
            
            # Check if this rule variant produces logs (clipping/sorting do not define logs)
            if stage in ["filter", "dup", "target"]:
                # Reconstruct the prefix string up to this point
                stages_prefix = current_string[:match.start()]
                
                # Recreate the exact {stages} value Snakemake generated for this step
                # Account for the temporary trailing underscore variant if present
                if match.group(3) == "_":
                    options_str = f"{options}_"
                else:
                    options_str = options
                    
                # Construct the paths to both log types defined in your rules
                collected_logs[f"{stage}{options_str}-{idx}"] = f"<sp>logs/{stages_prefix}-{stage}{options_str}-{idx}.stats.txt"

            # Move one step backward in the transformation grammar
            current_string = current_string[:match.start()]

    print("Collected logs:", collected_logs)
    return collected_logs


def construct_log_args(w, input):
    """
    Constructs the command-line arguments for the collect-stats tool
    based on the parsed log files from the final BAM filename.
    """
    args = []
    for key, log_path in input.items():
        if key in ["trims", "maps", "depths"]:
            continue  # Skip these keys as they are handled separately
        args.append(f"'{log_path}::{key}'")
    return ' '.join(args)


rule collect_pipeline_logs:
    """
    Aggregates all logs generated by the dynamic grammar path
    resolved from the 'final_bam_file' target configuration.
    """
    input:
        # Pass our recursive parser directly as an input function,
        unpack(parse_pipeline_logs),
        trims = expand_sp('<sp>logs/trimming_stat-{idx}.json'),
        maps = expand_sp('<sp>logs/mapped-{idx}.stats.txt'),
        depths = expand_sp('<sp>logs/mapped-final-{idx}.depths.txt.gz'),
    output:
        collected = '<sp>logs/stats.tsv'
    params:
        sample = get_sample,
        args = construct_log_args,
        trimmed = lambda wildcards, input: '--trimmed ' + ' --trimmed '.join(input.trims),
        mapped = lambda wildcards, input: '--mapped ' + ' --mapped '.join(input.maps),
        depthed = lambda wildcards, input: '--depth ' + ' --depth '.join(input.depths),

    threads: 1
    shell:
        'ngs-pl calculate-cpes-stats -o {output.collected} --mindepth {min_depth}'
        ' {params.trimmed} {params.mapped} {params.depthed} --sample {params.sample}'
        ' {params.args}'


# EOF
