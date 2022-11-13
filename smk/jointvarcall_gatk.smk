
rule jointvarcall_gatk:

    shell:
        "gatk GenotypeGVCFs -stand-call-conf 30 -R {refseq} --variant {gvcfs} -O joint-varcall.vcf"

