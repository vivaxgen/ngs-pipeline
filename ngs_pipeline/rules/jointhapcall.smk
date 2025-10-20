
# NOTE:
# to use this, use --target varcall_result

# we use Clair3 variant calling result and GATK GenotypeGVCFs to perform
# joint-variant calling

include: "jointvarcall_gatk.smk"

# directory notes
# SAMPLE_DIRS is [aboslute_path, absolute_path, ...]
# destdir is absolute_path to joint-variant callind result directory


rule haplotype_report:
    localrule: True
    input:
        #f"{destdir}/concatenated-bcsq.vcf.gz"
        f"{destdir}/haplotype-report.tsv"


def get_haplotype_files(w):
    # traversing on all source directories
    print(SAMPLE_DIRS)
    return [f'{a_dir}/haplotypes.tsv' for (a_dir, s) in zip(SAMPLE_DIRS, SAMPLES)]


rule prepare_haplotype:
    threads: 1
    input:
        bam = "{sample_dir}/maps/mapped-final.bam",
        vcf = f"{destdir}/concatenated-bcsq.vcf.gz",
    output:
        tsv = "{sample_dir}/haplotypes-nuc.tsv",
        vcf = "{sample_dir}/vcf/phased.vcf.gz",
    log:
        log1 = "{sample_dir}/logs/haplotype-report.json",
        log2 = "{sample_dir}/logs/construct-haplotypes.log"
    params:
        sample = subpath(subpath(output.tsv, ancestor=1), basename=True)
    shell:
        "ngs-pl construct-haplotypes"
        "  --sample {params.sample}"
        "  --min-match-length 2000 --outvcf {output.vcf}"
        "  --min-qual 10"
        "  --cds-only"
        "  --min-depth-ratio 0.5"
        "  --vcf-file {input.vcf}"
        "  --outlog {log.log1}"
        "  --outfile {output.tsv}"
        "  {input.bam}"
        "  2> {log.log2}"


rule reannotate_phased_vcf:
    threads: 1
    input:
        vcf = "{a_dir}/vcf/phased.vcf.gz",
    output:
        vcf = "{a_dir}/vcf/phased-bcsq-nonlocal.vcf.gz",
        tbi = "{a_dir}/vcf/phased-bcsq-nonlocal.vcf.gz.tbi"
    shell:
        "bcftools csq -o {output.vcf} -f {refseq} --gff {gff_file} {input.vcf}"
        " && bcftools index -t {output.vcf}"


rule refine_haplotype:
    threads: 1
    input:
        tsv = "{a_dir}/haplotypes-nuc.tsv",
        vcf = "{a_dir}/vcf/phased-bcsq-nonlocal.vcf.gz",
    output:
        tsv = "{a_dir}/haplotypes.tsv",
    shell:
        "ngs-pl recode-haplotypes"
        "  --phased-vcf {input.vcf}"
        "  --outfile {output.tsv}"
        "  {input.tsv}"


rule gather_haplotypes:
    threads: 1
    input:
        get_haplotype_files,
    output:
        f"{destdir}/haplotype-report.tsv"
    shell:
        "ngs-pl gather-haplotype"
        "  --outfile {output}"
        "  {input}"

# EOF
