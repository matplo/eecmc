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

if [ "$1" == "help" ] || [ -z "$1" ]; then
  echo "Usage: $0 <output_dir> [clean|jet_pt] |nev] [vincia|dire]"
  exit 0
fi

OUTPUT_DIR=$1
if [ -z ${OUTPUT_DIR} ]; then
  echo "Usage: $0 <output_dir>"
  exit 1
fi

if [ -d ${OUTPUT_DIR} ]; then
  if [ $2 == "clean" ]; then
    rm -rfv ${OUTPUT_DIR}/*
    exit 0
  fi
fi

JET_PT=20
if [ ! -z $2 ]; then
  JET_PT=$2
fi

NEV=100
if [ ! -z $3 ]; then
  NEV=$3
fi

PSHOWER=""
if [ ! -z $4 ]; then
  if [ "$4" == "vincia" ] || [ "$4" == "dire" ]; then
    PSHOWER="--py-$4"
    OUTPUT_DIR=${OUTPUT_DIR}/${4}
  else
    OUTPUT_DIR=${OUTPUT_DIR}/std
  fi
else
  OUTPUT_DIR=${OUTPUT_DIR}/std
fi

OUTPUT_DIR=${OUTPUT_DIR}/${JET_PT}

mkdir -p ${OUTPUT_DIR}
if [ ! -d ${OUTPUT_DIR} ]; then
  echo "Error: ${OUTPUT_DIR} does not exist"
  exit 1
fi

echo "SYS_YASP_DIR: ${SYS_YASP_DIR}"
echo "OUTPUT_DIR: ${OUTPUT_DIR}"
OUTPUT_FILE=${OUTPUT_DIR}/pythia8_jetreco_eec_0.root
JOB_SCRIPT=${OUTPUT_DIR}/_eec_pythia8_job_0.sh
COUNTER=200
while [ -e "${JOB_SCRIPT}" ]; do
  COUNTER=$((COUNTER+1))
  OUTPUT_FILE="${OUTPUT_DIR}/pythia8_jetreco_eec_${COUNTER}.root"
  JOB_SCRIPT=${OUTPUT_DIR}/_eec_pythia8_job_${COUNTER}.sh
done

# this is something to execute...
skeleton_script=${THIS_DIR}/slurm_eec_pythia8.sh
TMP_OUTPUT_DIR=/scratch/u/${OUTPUT_DIR}
TMP_OUTPUT_FILE=$(basename ${OUTPUT_FILE})
setup_cmnd="mkdir -p ${TMP_OUTPUT_DIR} && cd ${TMP_OUTPUT_DIR}"
cmnd="${THIS_DIR}/pythia8_jetreco_eec.py --py-pthatmin $JET_PT --nev ${NEV} --output ${TMP_OUTPUT_FILE} --py-hardQCD --py-seed -1 ${PSHOWER}"
cp_cmnd="cp -v ${TMP_OUTPUT_FILE} ${OUTPUT_FILE}"
yasprepl -f ${skeleton_script} -o ${JOB_SCRIPT} --define SYS_YASP_DIR=${SYS_YASP_DIR} OUTPUT_DIR=${OUTPUT_DIR} CMND_TO_RUN="${cmnd}" COPY_CMND="${cp_cmnd}" SETUP_CMND="${setup_cmnd}"

echo "Job script: ${JOB_SCRIPT}"
# cat ${JOB_SCRIPT}
sbatch ${JOB_SCRIPT}
echo "Job submitted"
echo "Output: ${OUTPUT_FILE}"
echo "Done"
exit 0
# End of file