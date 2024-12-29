#!/usr/bin/bash

# installation script for vivaxgen ngs-pipeline [https://github.com/vivaxgen/ngs-pipeline]

# optional variable:
# - BASEDIR
# - OMIT

set -eu

# run the base.sh
# Detect the shell from which the script was called
parent=$(ps -o comm $PPID |tail -1)
parent=${parent#-}  # remove the leading dash that login shells have
case "$parent" in
  # shells supported by `micromamba shell init`
  bash|fish|xonsh|zsh)
    shell=$parent
    ;;
  *)
    # use the login shell (basename of $SHELL) as a fallback
    shell=${SHELL##*/}
    ;;
esac

# Parsing arguments
if [ -t 0 ] && [ -z "${BASEDIR:-}" ]; then
  printf "Pipeline base directory? [./vvg-ngspl] "
  read BASEDIR
fi

# default value
BASEDIR="${BASEDIR:-./vvg-ngspl}"

uMAMBA_ENVNAME="${uMAMBA_ENVNAME:-ngs-pl}"
PYVER="${PYVER:-3.12}"
SNAKEMAKEVER="${SNAKEMAKEVER:-9}"
source <(curl -L https://raw.githubusercontent.com/vivaxgen/vvg-box/main/install.sh)

echo "Cloning vivaxGEN ngs-pipeline"
git clone --depth 1 https://github.com/vivaxgen/ngs-pipeline.git ${ENVS_DIR}/ngs-pipeline
ln -sr ${ENVS_DIR}/ngs-pipeline/etc/bashrc.d/10-ngs-pipeline ${BASHRC_DIR}/
ln -sr ${ENVS_DIR}/ngs-pipeline/etc/bashrc.d/95-prompt-history ${BASHRC_DIR}/

source ${ENVS_DIR}/ngs_pipeline/etc/inst-scripts/inst-deps.sh

echo "ngs-pipeline" >> ${ETC_DIR}/installed-repo.txt

echo
echo "vivaxGEN ngs-pipeline has been successfully installed. "
echo "Please read the docs for further setup."
echo "The base installation directory (VVG_BASEDIR) is:"
echo
echo `realpath ${BASEDIR}`
echo
echo "To activate the basic NGS-Pipeline environment (eg. for installing"
echo "or setting up base enviroment directory), execute the command:"
echo
echo `realpath ${BASEDIR}`/bin/activate
echo

# EOF
