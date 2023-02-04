# main activation script
# copy and modify this script to your NGSENV_BASEDIR directory
# parts to be check and modified are mentioned with: -- modify as necessary


_script="$(readlink -f ${BASH_SOURCE[0]})"

## Delete last component from $_script ##
_mydir="$(dirname $_script)"
echo "${_mydir}"

## Set NGSENV_BASEDIR
export NGSENV_BASEDIR=${_mydir}

## set JOBCMD to be used, leave empty for single server system
## -- modify as ncessary
export JOBCMD='sbatch --cpus-per-task={threads} --job-name=smk-{rule}-{sample} --output=slurm/j-%j.out'

## prepare cache directory for samtools fastq cram-to-fastq conversion
## see here: https://www.htslib.org/workflow/cram.html
export REF_PATH=${NGSENV_BASEDIR}/cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s
export REF_CACHE=${NGSENV_BASEDIR}/cache/%2s/%2s/%s

## run other activation scripts here (such as conda etc) 
## -- modify as necessary
. ${NGSENV_BASEDIR}/../opt/etc/bashrc

## run ngs-pipeline activation script
## -- modify as necessary
. ${NGSENV_BASEDIR}/../ngs-pipeline/bin/activate.sh

## set prompt here
## -- modify as necessary
PS1="(pv-wgs) [\u@\h \W]\$ "

# EOF
