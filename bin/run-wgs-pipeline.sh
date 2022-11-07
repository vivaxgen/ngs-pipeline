#!/bin/sh

# this will run the pipeline using batch job submission system

_script="$(readlink -f ${BASH_SOURCE[0]})"

## Delete last component from $_script ##
_mydir="$(dirname $_script)"
echo "${_mydir}"

# JOBCMD='srun -N 1 -t 48:00:00'

parallel --eta -j 16 --workdir $PWD/{} "${JOBCMD} ${_mydir}/run_varcall.py" ::: `ls`

# EOF
