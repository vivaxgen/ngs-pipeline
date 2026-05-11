# stage-2

INST_SCRIPTS_DIR="${ENVS_DIR}/ngs-pipeline/etc/inst-scripts"

${VVGBIN}/link-resource-files.sh ${ENVS_DIR}/ngs-pipeline/etc/bashrc.d

if [[ -z ${VVG_MANIFEST_FILE:-} ]]; then
  echo -e "\e[32m>>> No manifest file provided, installing dependencies with inst-deps.sh\e[0m"
  source ${INST_SCRIPTS_DIR}/inst-deps.sh
fi

# EOF