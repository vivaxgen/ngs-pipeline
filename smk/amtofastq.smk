
# this smk will create paired fastq files from any alignment map files
# that are recognized by samtools (SAM/BAM/CRAM)

SOURCES = config['SOURCES']


def get_fastq_files():
    return [f"{sample}_R1.fastq.gz" for sample in SOURCES.keys()]


def get_source_filename(sample):
    return SOURCES[sample]


localrules: create_manifest_file


rule all:
    input:
        "manifest.tsv"


# convert an AM file to 
rule am_to_fastq:
    threads: 2
    input:
        lambda w: get_source_filename(w.sample)
    output:
        read1 = '{sample}_R1.fastq.gz',
        read2 = '{sample}_R2.fastq.gz'
    shell:
        "samtools collate -u -f -O {input} | "
        "samtools fastq -1 {output.read1} -2 {output.read2} -0 /dev/null -s /dev/null -n"


# check we have all fastq files from each bam and create a manifest file
rule create_manifest_file:
    threads: 2
    input:
        lambda w: get_fastq_files()
    output:
        "manifest.tsv"
    shell:
        # write sample codes and their corresponding fastq files to manifest file
        "touch {output}"

# EOF
