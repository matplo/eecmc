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

SYS_YASP_DIR=$(yaspenv.sh yasp -q feature yasp_dir 2>&1 | tail -n 1)
source ${SYS_YASP_DIR}/venvyasp/bin/activate

module use $SYS_YASP_DIR/software/modules
module load moduleName >/dev/null 2>&1
# module list

OUTPUT_DIR=$1
if [ -z ${OUTPUT_DIR} ]; then
  echo "Usage: $0 <output_dir>"
  exit 1
fi

skeleton_script=${THIS_DIR}/slurm_eec_pythia8.sh
job_script=${OUTPUT_DIR}/eec_pythia8.sh
yasprepl -f ${skeleton_script} -o ${job_script} \\
	--define SYS_YASP_DIR=$SYS_YASP_DIR \\
  --define OUTPUT_DIR=$OUTPUT_DIR

cat ${job_script}