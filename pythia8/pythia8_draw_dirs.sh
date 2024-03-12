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
module load yasp root >/dev/null 2>&1


if [ "$1" == "help" ] || [ -z "$1" ]; then
  echo "Usage: $0 <output_dir>"
  exit 0
fi

OUTPUT_DIR=$1
if [ -z ${OUTPUT_DIR} ]; then
  echo "Usage: $0 <output_dir>"
  exit 1
fi

output_types="h_add_types.txt"

files=$(find ${OUTPUT_DIR} -name "pythia8_jetreco_eec*.root")
for fname in ${files}
do
	dirname=$(dirname "${fname}")
	jetpt=$(basename ${dirname})
	pshower=$(dirname ${dirname})
	pshower=$(basename ${pshower})
	basename=$(basename "${fname}")
	output="${dirname}/h_${basename}"
	hadd_file="${pshower}_${jetpt}.txt"
	[ "${jetpt}" == "20" ] && hname="hJetE2C_norm_per_bin_ch_jetR4_2040_trk10_pythia"
	[ "${jetpt}" == "40" ] && hname="hJetE2C_norm_per_bin_ch_jetR4_4060_trk10_pythia"
	[ "${jetpt}" == "60" ] && hname="hJetE2C_norm_per_bin_ch_jetR4_6080_trk10_pythia"
	hfname="${THIS_DIR}/wenqing_mc_output.root"
	if [ -f "${fname}" ]; then
		cmnd="${THIS_DIR}/draw_RL.py -i ${fname} --use-h "${hfname}:${hname}" --ncorrel 2 -o ${output}" 
		if [ ! -e ${output} ]; then
			echo "input file: ${fname} output file: ${output} jetpt: ${jetpt} hname: ${hname} pshower: ${pshower}"
			$cmnd
		fi
		# echo ${output} | tee -a ${hadd_file}
		echo ${fname} | tee -a ${hadd_file}
		echo ${hadd_file} | tee -a ${output_types}
	fi
done

h_add_types=$(cat ${output_types} | uniq)
for hadd_file in ${h_add_types}
do
	_rootfile="${hadd_file%.*}.root"
	echo "hadd -f ${_rootfile} @${hadd_file}"
	hadd -f ${_rootfile} @${hadd_file}
	jet_pt=$(echo ${_rootfile} | cut -d_ -f 2 | cut -d. -f 1)
	[ "${jetpt}" == "20" ] && hname="hJetE2C_norm_per_bin_ch_jetR4_2040_trk10_pythia"
	[ "${jetpt}" == "40" ] && hname="hJetE2C_norm_per_bin_ch_jetR4_4060_trk10_pythia"
	[ "${jetpt}" == "60" ] && hname="hJetE2C_norm_per_bin_ch_jetR4_6080_trk10_pythia"
	hfname="${THIS_DIR}/wenqing_mc_output.root"
	output="h_${_rootfile}"
	${THIS_DIR}/draw_RL.py -i ${_rootfile} --use-h "${hfname}:${hname}" --ncorrel 2 -o ${output}
done
