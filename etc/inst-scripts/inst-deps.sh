# install dependencies for ngs-pipeline

OMIT="${OMIT:-}"

if ! [[ "$OMIT" =~ GATK ]]; then
  echo "Installing the latest GATK"
  retry 5 micromamba -y install "GATK4>=4.6,<5" -c conda-forge -c bioconda
fi

echo "Installing other dependencies with micromamba"
retry 5 micromamba -y install -n ${uMAMBA_ENVNAME} -f ${ENVS_DIR}/ngs-pipeline/etc/inst-scripts/env.yaml
