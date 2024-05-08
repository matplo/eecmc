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
THISD=$(thisdir)

cd ${THISD}

savedir=${PWD}
nev=1000

if [ "$1" == '-h' ]; then
	echo "Usage: $0 [jetpt] [nev] [njobs] [lund]"
	exit 0
fi

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

echo "[i] using ${SYS_YASP_DIR}"
source ${SYS_YASP_DIR}/venvyasp/bin/activate
module use ${SYS_YASP_DIR}/software/modules
module avail
module load yasp
module load mysherpa

datfile_default=${THISD}/charmRun.dat
datfile=${datfile_default}
if [ "$4" == "lund" ]; then
	datfile="${THISD}/charmRunDISLundTune.dat"
	dirmod="_lund"
fi

slurm_script=${THISD}/sherpa_slurm_job.sh

dname_base="/rstorage/ploskon/eec_sherpa_charm"
mkdir -pv ${dname_base}

jetpt=15
if [ ! -z $1 ]; then
	jetpt=$1
fi

if [ ! -z $2 ]; then
	nev=$2
fi

njobs=1
if [ ! -z $3 ]; then
	njobs=$3
fi

echo "[i] jetpt=${jetpt} nev=${nev} njobs=${njobs} datfile=${datfile}"
SDATE=$(date +"%Y%m%d%H%M")
jobs_to_run_shell_submit=${THISD}/submit_jobs_to_run_${SDATE}.sh
jobs_to_run_shell=${THISD}/jobs_to_run_${SDATE}.sh
for nj in `seq 1 ${njobs}`
do
	rseed=1000000
	dname="${dname_base}/charm_jetpt${jetpt}${dirmod}/${rseed}"
	while [ -d "${dname}" ]; do
		rseed=$((rseed-1))
		dname="${dname_base}/charm_jetpt${jetpt}${dirmod}/${rseed}"
	done
	dname="${dname_base}/charm_jetpt${jetpt}${dirmod}/${rseed}"
	if [ -d ${dname} ]; then
		echo "[i] directory exists - skipping - ${dname}"
	else
		echo "[i] creating directory - ${dname}"
		mkdir -p ${dname}
		cd ${dname}
		yasprepl -f ${datfile} \
					--define \
					jet_min_pt="${jetpt}.0" \
					jet_eta_max=0.5 \
					jet_R=0.4 \
					p_beam_energy=2510 \
					random_seed=${rseed} \
					-o ${dname}/Run.dat
		yasprepl -f ${slurm_script} \
					--define \
					number_of_events=${nev} \
					output_dir=${dname} \
					-o ${dname}/sherpa_slurm_job.sh
		chmod +x sherpa_slurm_job.sh
		chmod +x ${slurm_script}
		echo "${dname}/sherpa_slurm_job.sh" >> ${jobs_to_run_shell}
		chmod +x ${jobs_to_run_shell}
		ln -sf ${jobs_to_run_shell} ${THISD}/latest_jobs_to_run.sh
		echo "[i] run latest_jobs_to_run.sh to run sequentially"

		echo "sbatch ${dname}/sherpa_slurm_job.sh" >> ${jobs_to_run_shell_submit}
		chmod +x ${jobs_to_run_shell_submit}
		ln -sf ${jobs_to_run_shell_submit} ${THISD}/submit_latest_jobs_to_run.sh
		echo "[i] run submit_latest_jobs_to_run.sh to submit jobs"

		analysis_config=${THISD}/analysis_config_D_jetpt${jetpt}${dirmod}.yaml
		cp -v ${analysis_config} ./analysis_config.yaml
		cp -v ${THISD}/heec.py ${dname}/
		cp -v ${THISD}/hepmc_D_analysis.py ${dname}/
		cp -v ${THISD}/hepmc_count_events.py ${dname}/
		cd ${savedir}
	fi
done

cd ${savedir}
