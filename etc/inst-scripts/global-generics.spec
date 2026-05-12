# generics
bcftools>=1.23 samtools>=1.23 vcftools>=0.1.16  --environment generic-tools
fastp>=1.3
bwa-mem2>=2.3 bwa>=0.7.19 minimap2>=2.29 bowtie2>=2.5.4 --environment mapper
freebayes>=1.3.10
snpeff>5.3 fastqc>=0.12.1 ${JAVA_BASED} --environment java-based
sambamba>=1.0.1
${OPTIONAL_PACKAGES}


# EOF
