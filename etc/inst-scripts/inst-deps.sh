# install dependencies for ngs-pipeline

echo -e "\e[32m>>>> Setting pixi default channels to conda-forge and bioconda\e[0m"
pixi config set default-channels '["conda-forge", "bioconda"]' --global
pixi config set default-channels '["conda-forge", "bioconda"]'

if ! defined_and_contains_any VVG_EXCLUDE gatk4; then
  echo -e "\e[32m>>>> Installing the latest gatk4\e[0m"
  pixi-global-install ${INST_SCRIPTS_DIR}/global-gatk4.spec
fi

echo -e "\e[32m>>>> Installing global generic dependencies\e[0m"
pixi-global-install ${INST_SCRIPTS_DIR}/global-generics.spec

# if PIXI_ENVIRONMENT_PLATFORMS is defined and there is exist filename
# ${PIXI_ENVIRONMENT_PLATFORMS}.spec, then install the dependencies in that file
if [[ -n "${PIXI_ENVIRONMENT_PLATFORMS:-}" ]] && [ -f ${INST_SCRIPTS_DIR}/global-${PIXI_ENVIRONMENT_PLATFORMS}.spec ]; then
  echo -e "\e[32m>>>> Installing global dependencies for platform ${PIXI_ENVIRONMENT_PLATFORMS}\e[0m"
  pixi-global-install ${INST_SCRIPTS_DIR}/global-${PIXI_ENVIRONMENT_PLATFORMS}.spec
fi

# installing workspace-specific dependencies
echo -e "\e[32m>>>> Installing workspace-specific dependencies\e[0m"
pixi-add ${INST_SCRIPTS_DIR}/workspace-generics.spec

# EOF
