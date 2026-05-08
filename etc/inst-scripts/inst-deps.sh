# install dependencies for ngs-pipeline

OMIT="${OMIT:-}"

if ! [[ "$OMIT" =~ GATK ]]; then
  echo "Installing the latest GATK"
  retry 5 micromamba -vvv -y install "GATK4>=4.6,<5" -c conda-forge -c bioconda -c defaults
fi

echo "Installing other dependencies with micromamba"
retry 5 micromamba -vvv -y install -n ${uMAMBA_ENVNAME} -f ${ENVS_DIR}/ngs-pipeline/etc/inst-scripts/env.yaml

echo "Installing Clair3 as ngs-pl-clair3 sub-environment"
retry 5 micromamba -vvv create -y -n ${uMAMBA_ENVNAME}-clair3 clair3 python=3.9.0 -c conda-forge -c bioconda -c defaults


