
# prepare necessary parameters

import os

ngs_pipeline_basedir = os.environ['NGS_PIPELINE_BASE']
refmap = config.get('refmap_file', None)
refseq = config['refseq_file']


# check available read files

wildcard_constraints:
    idx = '\\d+'

IDXS, = glob_wildcards('reads/raw-{idx}_R1.fastq.gz')

# final output of this workflow
rule all:
    input:
        "maps/mapped.bam"

include: config.get('reads_trimmer_wf', 'trimmer_cutadapt.smk')
include: config.get('reads_mapper_wf', 'mapper_minimap2.smk')

rule map_merging:
    input:
        expand('maps/mapped-{idx}.bam', idx=IDXS)
    output:
        'maps/mapped.bam'
    run:
        if len(IDXS) > 1:
            shell('samtools merge -@8 {output} {input}')
        else:
            shell('ln ln -s {input} {output}')

# EOF
