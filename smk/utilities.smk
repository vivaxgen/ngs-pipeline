
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
        fasta = "{pfx}/{fn}.fasta.fai"
    shell:
        "samtools faidx {input.fasta}"


rule index_ont_mmi:
    threads: 1
    input:
        fasta = "{pfx}/{fn}.fasta"
    output:
        index = "{pfx}/{fn}.fasta.ont.mmi"
    shell:
        "minimap2 -x map-ont -d {output.index} {input.fasta}"


rule bunzip2:
    threads: 1
    input:
        path = "{fn}.bz2"
    output:
        path = "{fn}"
    shell:
        "bunzip2 -c {input.path} > {output.path}"


# EOF
