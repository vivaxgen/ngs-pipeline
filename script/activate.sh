
## source this first to set up proper environment

_script="$(readlink -f ${BASH_SOURCE[0]})"

## Delete last component from $_script ##
_mydir="$(dirname $_script)"

export NGS_PIPELINE_BASE="$(dirname $_mydir)"

