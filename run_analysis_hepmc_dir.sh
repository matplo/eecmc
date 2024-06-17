#!/bin/bash

input_dir=$1
if [ -z ${input_dir} ]; then
	echo "Usage: $0 <input_dir>"
	exit 1
fi
if [ -d ${input_dir} ]; then
	# find all hepmc files
	echo "Finding hepmc files in ${input_dir}"
	hepmc_files=$(find ${input_dir} -name "*.hepmc")
else
	if [ -f ${input_dir} ]; then
		hepmc_files=$(cat ${input_dir})
	else
		echo "Usage: $0 <input_dir>"
		exit 1
	fi
fi

if [ ${#hepmc_files[@]} -lt 1 ]; then
	echo "Error: hepmc_files not found"
	exit 1
fi

echo "List has ${#hepmc_files[@]} hepmc files"

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
export -f thisdir

THIS_DIR=$(thisdir)

function run_on_a_file()
{
	input_file=$1
	EECMC_DIR=/software/users/ploskon/eecmc
	if [ -z ${input_file} ]; then
		echo "Error: input_file not found"
		return 1
	fi
	echo "Running on ${input_file}"
	# setup jet pt cuts
	if [[ ${input_file} == *"jetpt10"* ]]; then
		jet_min_pt=10
		jet_max_pt=15
		if [[ ${input_file} == *"5020"* ]]; then
			jet_max_pt=20
		fi
	else
		jet_min_pt=15
		jet_max_pt=30
	fi
	charm_setting=""
	# check if charm or inclusive
	if [[ ${input_file} == *"charm"* ]]; then
		echo "Running on charm"
		charm_setting="--D0mode 1 --D0-pt-min 5 --D0-pt-max ${jet_max_pt}"
	else
		echo "Running on inclusive"
	fi
	jet_setting="--jet-pt-min ${jet_min_pt} --jet-pt-max ${jet_max_pt}"
	echo "Jet setting: ${jet_setting}"
	echo "Charm setting: ${charm_setting}"

	# setup the environment
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
	module load bundle/hepbase
	module load heppyy/current

	if [ ! -e $EECMC_DIR/eecmc.module ]; then
		echo "Error: $EECMC_DIR/eecmc.module not found"
		exit 1
	fi
	module use $EECMC_DIR
	module load eecmc.module

	output_file=$(dirname ${input_file})/$(basename ${input_file} .hepmc)_ana_eec.root
	script_file=$(dirname ${input_file})/$(basename ${input_file} .hepmc)_ana_eec.sh
	log_file=$(dirname ${input_file})/$(basename ${input_file} .hepmc)_ana_eec.log
	nev=10000
	echo "simple_jet_rec.py --enable-eec --nev ${nev} --output ${output_file} --hepmc 3 --input ${input_file} ${jet_setting} ${charm_setting}" | tee ${script_file}
	chmod +x ${script_file}
	# ${script_file} 2>&1 | tee ${log_file}
	simple_jet_rec.py --enable-eec --nev ${nev} --output ${output_file} --hepmc 3 --input ${input_file} ${jet_setting} ${charm_setting} 2>&1 | tee ${log_file}
}
export -f run_on_a_file

# use gnu parallel to run on multiple files with max 5 jobs
# run gnuparallel with 5 jobs max on each file run_on_a_file with status how many to go
# run_on_a_file will run the analysis on the file
# gnuparallel will run the analysis on the file

parallel -j 1 --bar --eta run_on_a_file ::: ${hepmc_files}

# # run on first few files
# for file in ${hepmc_files}; do
# 	run_on_a_file ${file}
# 	break
# done
