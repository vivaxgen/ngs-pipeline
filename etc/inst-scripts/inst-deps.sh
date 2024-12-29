# install dependencies for ngs-pipeline

OMIT="${OMIT:-}"

if ! [[ "$OMIT" =~ GATK ]]; then
  echo "Installing the latest GATK"
  retry 5 micromamba -y install GATK4 -c conda-forge -c bioconda
fi

echo "Installing latest htslib tools"
micromamba -y install "bcftools>=1.18" "samtools>=1.18" -c conda-forge -c bioconda -c defaults

echo "Installing latest fastp"
micromamba -y install fastp -c conda-forge -c bioconda

echo "Installing latest bwa-mem2"
micromamba -y install bwa-mem2 -c conda-forge -c bioconda

echo "Installing latest bwa"
micromamba -y install bwa -c conda-forge -c bioconda

echo "Installing minimap2"
micromamba -y install minimap2 -c conda-forge -c bioconda -c defaults

echo "Installing bowtie2"
micromamba -y install bowtie2 -c conda-forge -c bioconda -c defaults

echo "Installing freebayes"
micromamba -y install freebayes=1.3.6 -c conda-forge -c bioconda -c defaults

echo "Installing vcftools"
micromamba -y install vcftools -c conda-forge -c bioconda

echo "Installing chopper"
micromamba -y install chopper -c conda-forge -c bioconda -c defaults

echo "Installing snpEff"
micromamba -y install snpeff -c conda-forge -c bioconda -c defaults

echo "Installing sambamba"
micromamba -y install sambamba -c conda-forge -c bioconda -c defaults

echo "Installing fastqc"
micromamba -y install fastqc -c conda-forge -c bioconda -c defaults

echo "Installing required Python modules"
pip3 install "snakemake<${SNAKEMAKEVER}"
pip3 install snakemake-executor-plugin-cluster-generic
pip3 install cyvcf2
pip3 install pysam
pip3 install pandas
pip3 install Pillow
pip3 install IPython
pip3 install matplotlib
pip3 install multiqc
pip3 install cutadapt

# pip3 install pycairo
# we use conda pycairo since pip pycairo does not have complete binary
# distribution
micromamba -y install pycairo -c conda-forge

pip3 install NanoPlot
pip3 install argcomplete
pip3 install openpyxl

# EOF