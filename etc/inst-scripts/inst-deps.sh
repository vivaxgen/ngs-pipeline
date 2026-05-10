# install dependencies for ngs-pipeline

INST_SCRIPTS_DIR="${ENVS_DIR}/ngs-pipeline/etc/inst-scripts"

${VVGBIN}/link-resource-files.sh ${ENVS_DIR}/ngs-pipeline/etc/bashrc.d

echo "Setting pixi default channels to conda-forge and bioconda"
pixi config set default-channels '["conda-forge", "bioconda"]' --global
pixi config set default-channels '["conda-forge", "bioconda"]'

if ! defined_and_contains_any VVG_EXCLUDE gatk4; then
  echo "Installing the latest gatk4"
  pixi-global-install ${INST_SCRIPTS_DIR}/global-gatk4.spec
fi

echo "Installing global generic dependencies"
pixi-global-install ${INST_SCRIPTS_DIR}/global-generics.spec

# if PIXI_ENVIRONMENT_PLATFORMS is defined and there is exist filename
# ${PIXI_ENVIRONMENT_PLATFORMS}.spec, then install the dependencies in that file
if [[ -n "${PIXI_ENVIRONMENT_PLATFORMS:-}" ]] && [ -f ${INST_SCRIPTS_DIR}/global-${PIXI_ENVIRONMENT_PLATFORMS}.spec ]; then
  echo "Installing global dependencies for platform ${PIXI_ENVIRONMENT_PLATFORMS}"
  pixi-global-install ${INST_SCRIPTS_DIR}/global-${PIXI_ENVIRONMENT_PLATFORMS}.spec
fi
# installing workspace-specific dependencies
echo "Installing workspace-specific dependencies"
pixi-add ${INST_SCRIPTS_DIR}/workspace-generics.spec

# EOF
