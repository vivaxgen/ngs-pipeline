# install dependencies for ngs-pipeline

echo -e "\e[32m>>>> Setting pixi default channels to conda-forge and bioconda\e[0m"
pixi config set default-channels '["conda-forge", "bioconda"]' --global
pixi config set default-channels '["conda-forge", "bioconda"]'

JAVA_BASED=""
OPTIONAL_PACKAGES=""
if ! defined_and_contains_any VVG_EXCLUDE gatk4; then
  echo -e "\e[32m>>>> Will install the latest gatk4\e[0m"
  JAVA_BASED="${JAVA_BASED} gatk4=4.6"
  #pixi-global-install ${INST_SCRIPTS_DIR}/global-gatk4.spec
fi


LINUX64_ONLY=""
if ! defined_and_contains_any VVG_EXCLUDE ONT_TOOLS; then
  echo -e "\e[32m>>>> Will install the latest ONT data tools\e[0m"
  OPTIONAL_PACKAGES="${OPTIONAL_PACKAGES} fastplong>=0.4.1 clair3>=2.0"
  LINUX64_ONLY="glnexus>=1.4"
fi

if [[ -n "${PIXI_ENVIRONMENT_PLATFORMS:-}"

echo -e "\e[32m>>>> Installing global generic dependencies\e[0m"
pixi-global-install ${INST_SCRIPTS_DIR}/global-generics.spec

# if PIXI_ENVIRONMENT_PLATFORMS is defined and there is exist filename
# ${PIXI_ENVIRONMENT_PLATFORMS}.spec, then install the dependencies in that file
if [[ -n "${PIXI_ENVIRONMENT_PLATFORMS:-}" ]] && [ -f ${INST_SCRIPTS_DIR}/global-${PIXI_ENVIRONMENT_PLATFORMS}.spec ]; then
  if [[ "${PIXI_ENVIRONMENT_PLATFORMS}" == "linux-64" ]]; then
    OPTIONAL_PACKAGES="${OPTIONAL_PACKAGES:+$OPTIONAL_PACKAGES }$LINUX64_ONLY"
  fi
  echo -e "\e[32m>>>> Installing global dependencies for platform ${PIXI_ENVIRONMENT_PLATFORMS}\e[0m"
  pixi-global-install ${INST_SCRIPTS_DIR}/global-${PIXI_ENVIRONMENT_PLATFORMS}.spec
fi

# installing workspace-specific dependencies
echo -e "\e[32m>>>> Installing workspace-specific dependencies\e[0m"
pixi-add ${INST_SCRIPTS_DIR}/workspace-generics.spec

# EOF
