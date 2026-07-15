

indir = config.get("indir").removesuffix("/")
outdir = config.get("outdir").removesuffix("/")
if outdir[0] != '/':
    outdir = f"{indir}/{outdir}"

print(f"Outdir is: {outdir}")

fastqc_files = [f"{outdir}/{fn}_fastqc.zip" for fn in glob_wildcards(f"{indir}/{{fn}}.fastq.gz")[0]]
print(fastqc_files)


rule all:
    input:
        f"{outdir}/multiqc_report.html"


rule run_fastqc:
    threads: 2
    input:
        fastq = f"{indir}/{{fn}}.fastq.gz"
    output:
        fastqc = f"{outdir}/{{fn}}_fastqc.zip"
    shell:
        "fastqc {input.fastq} -o {outdir}"


rule run_multiqc:
    threads: 1
    input:
        fastqc_files,
    output:
        multiqc = f"{outdir}/multiqc_report.html"
    shell:
        "multiqc -o {outdir} {outdir}"


# EOF