
# helper rules to be used across mappers

rule mapped_bam_link:
    # this rule is used if the mapped bam file needs to contain sample name
    # eg. to make the file unique across multiple samples
    localrule: True
    input:
        bam = "<sp>maps/mapped-{idx}.bam",
    output:
        bam = temp_unless(get_mapped_bam_file(), keep_paired_bam),
    shell:
        "ln -f {input.bam} {output.bam}"


# EOF
