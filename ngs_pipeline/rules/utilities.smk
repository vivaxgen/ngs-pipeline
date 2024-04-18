
# utilities

rule index_csi:
    threads: 1
    input:
        "{pfx}/{reg}.vcf.gz"
    output:
        "{pfx}/{reg}.vcf.gz.csi"
    shell:
        "bcftools index --csi {input}"


rule index_tbi:
    threads: 1
    input:
        "{pfx}/{reg}.vcf.gz"
    output:
        "{pfx}/{reg}.vcf.gz.tbi"
    shell:
        "bcftools index --tbi {input}"


rule index_bai:
    threads: 1
    input:
        "{pfx}/{reg}.bam"
    output:
        "{pfx}/{reg}.bam.bai"
    shell:
        "samtools index {input}"


rule index_fai:
    threads: 1
    input:
        fasta = "{pfx}/{fn}.fasta"
    output:
        fai = "{pfx}/{fn}.fasta.fai"
    shell:
        "samtools faidx {input.fasta}"


rule index_dict:
    threads: 1
    input:
        fasta = "{pfx}/{fn}.fasta"
    output:
        dict = "{pfx}/{fn}.dict"
    shell:
        "samtools dict -o {output.dict} {input.fasta}"


rule index_ont_mmi:
    threads: 1
    input:
        fasta = "{pfx}/{fn}.fasta"
    output:
        index = "{pfx}/{fn}.fasta.ont.mmi"
    resources:
        mem_mb = config.get('index_mem_gb', 16) * 1024
    shell:
        "minimap2 -x map-ont -d {output.index} {input.fasta}"


rule index_bwamem2:
    threads: 1
    input:
        fasta = "{pfx}/{fn}.fasta"
    output:
        index = "{pfx}/{fn}.fasta.bwt.2bit.64"
    resources:
        mem_mb = config.get('index_mem_gb', 16) * 1024
    shell:
        "bwa-mem2 index {input.fasta}"


rule index_bwa:
    threads: 1
    input:
        fasta = "{pfx}/{fn}.fasta"
    output:
        index = "{pfx}/{fn}.fasta.bwt"
    resources:
        mem_mb = config.get('index_mem_mb', 16384)
    shell:
        "bwa index {input.fasta}"


rule bunzip2:
    threads: 1
    input:
        path = "{fn}.fasta.bz2"
    output:
        path = "{fn}.fasta"
    shell:
        "bunzip2 -c {input.path} > {output.path}"


# EOF
