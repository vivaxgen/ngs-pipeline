#!/usr/bin/bash

# installation script for vivaxgen ngs-pipeline [https://github.com/vivaxgen/ngs-pipeline]

# optional variable:
# - VVG_BASEDIR
# - PIXI_ENVNAME
# - VVG_EXCLUDE
# - VVG_INCLUDE
# - VVG_NGSPL_REPOURL

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
if [ -t 0 ] && [ -z "${VVG_BASEDIR:-}" ]; then
  printf "Pipeline base directory? [./vvg-ngspl] "
  read VVG_BASEDIR
fi

# default value
VVG_BASEDIR="${VVG_BASEDIR:-./vvg-ngspl}"

PIXI_ENVNAME="${PIXI_ENVNAME:-ngs-pl}"

# install pixi-based vvv-box
echo ">> Installing pixi-based vvg-box"
source <(curl -L https://raw.githubusercontent.com/vivaxgen/vvg-box/main/install.sh)

echo ">> Cloning vivaxGEN ngs-pipeline repository"
# add --branch dev for dev
git clone --depth 1 ${VVG_NGSPL_REPOURL:-https://github.com/vivaxgen/ngs-pipeline.git} ${ENVS_DIR}/ngs-pipeline  

# source the 2nd stage installation script for dependencies
echo ">> Executing stage 2 of installation"
source ${ENVS_DIR}/ngs-pipeline/etc/inst-scripts/inst-deps.sh

echo ">> Adding ngs-pipeline to installed repositories"
echo "ngs-pipeline" >> ${ETC_DIR}/installed-repo.txt

echo
echo "vivaxGEN ngs-pipeline has been successfully installed. "
echo "Please read the docs for further setup."
echo "The base installation directory (VVG_BASEDIR) is:"
echo
echo `realpath ${VVG_BASEDIR}`
echo
echo "To activate the basic NGS-Pipeline environment (eg. for installing"
echo "or setting up base enviroment directory), execute the command:"
echo
echo `realpath ${VVG_BASEDIR}`/bin/activate
echo

# EOF
