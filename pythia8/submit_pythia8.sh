#!/bin/bash

function thisdir()
{
        SOURCE="${BASH_SOURCE[0]}"
        while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
          DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
          SOURCE="$(readlink "$SOURCE")"
          [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
        done
        DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
        echo ${DIR}
}
THIS_DIR=$(thisdir)
yaspenv_shell=$(which yaspenv.sh)
if [ -z ${yaspenv_shell} ]; then
  echo "Error: yaspenv.sh not found"
  exit 1
fi
SYS_YASP_DIR=$(${yaspenv_shell} yasp -q feature yasp_dir 2>&1 | tail -n 1)
if [ -z ${SYS_YASP_DIR} ]; then
  echo "Error: SYS_YASP_DIR not found"
  exit 1
fi
source ${SYS_YASP_DIR}/venvyasp/bin/activate

module use $SYS_YASP_DIR/software/modules
module load yasp >/dev/null 2>&1
# module list

OUTPUT_DIR=$1
if [ -z ${OUTPUT_DIR} ]; then
  echo "Usage: $0 <output_dir>"
  exit 1
fi

mkdir -p ${OUTPUT_DIR}
if [ ! -d ${OUTPUT_DIR} ]; then
  echo "Error: ${OUTPUT_DIR} does not exist"
  exit 1
fi

echo "SYS_YASP_DIR: ${SYS_YASP_DIR}"
echo "OUTPUT_DIR: ${OUTPUT_DIR}"

skeleton_script=${THIS_DIR}/slurm_eec_pythia8.sh
job_script=${OUTPUT_DIR}/_eec_pythia8_job.sh
yasprepl -f ${skeleton_script} -o ${job_script} --define SYS_YASP_DIR=${SYS_YASP_DIR} OUTPUT_DIR=${OUTPUT_DIR}

echo "Job script: ${job_script}"
cat ${job_script}
