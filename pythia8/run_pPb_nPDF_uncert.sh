#!/bin/bash

# replace the nPDF setup in the config file with another number from the set
function replace_nPDF()
{
		local cmnd_file=$1
		local nPDFname=$2
		local nPDF=$3
		local output_dir=$4
		_cmnd_file=$(basename ${cmnd_file})
		local output_file=${output_dir}/${_cmnd_file%.cmnd}_nPDF_${nPDFname}_${nPDF}.cmnd

		if [ ! -f $cmnd_file ]; then
				echo "File $cmnd_file does not exist."
				exit 1
		fi

		# sed -e "s|.*EPPS21nlo_CT18Anlo_Pb208.*|${nPDFname}:${nPDF}|g" $cmnd_file > $output_file
		# replace line with "PDF:pSetB = ..." with "PDF:pSetB = LHAPDF6:${nPDFname}:${nPDF}"
		sed -e "s|.*PDF:pSetB.*|PDF:pSetB = LHAPDF6:${nPDFname}/${nPDF}|g" $cmnd_file > $output_file
		echo $output_file
}
export -f replace_nPDF


# run Pythia with the given command file
run_pythia() {
    local cmnd_file=$1
    local nev=$2
    local charged=$3
    # local output_file="simple_eec_${cmnd_file%.cmnd}_output.root"
    local output_file="${cmnd_file%.cmnd}.root"
		if [ ! -z $4 ]; then
				output_file=$4
		fi
    
    if [ ! -f $cmnd_file ]; then
        echo "File $cmnd_file does not exist."
        exit 1
    fi

		# for reproducibility ...
		local pythia_seed=123456

		if [ -f ${output_file} ]; then
				echo "File ${output_file} already exists. Skipping."
				return
		fi
    ./pythia8_simple_eec.py --py-cmnd ${cmnd_file} --nev ${nev} --jet-pt-min 20. ${charged} --py-pthatmin 20. --py-seed ${pythia_seed} -o ${output_file} 2>&1 | tee ${output_file%.root}.log
}
export -f run_pythia

. ./util.sh

help_requested=$(get_opt "help" $@)
if [ "x${help_requested}" == "xyes" ]; then
		echo_info "Usage: $0 [options]"
		echo_info "Options:"
		echo_info "  --nPDFset=<nPDFset> : nPDF set to use (default: EPPS21nlo_CT18Anlo_Pb208)"
		echo_info "  --cmnd=<cmnd_file> : command file to use"
		echo_info "  --output_prefix=<output_prefix> : prefix for output files (default: \$PWD/nPDFgen)"
		echo_info "  --output_dir=<output_dir> : output directory (default: <output_prefix>/<cmnd_file>-<nPDFset>)"
		echo_info "  --nev=<nev> : number of events to generate (default: 100000)"
		echo_info "  --nproc=<nproc> : number of processes to run in parallel (default: 5)"
		exit 0
fi

nev=$(get_opt "nev" $@)
if [ -z ${nev} ]; then
		nev=10000
		echo_note "Number of events not specified. Using default: ${nev}"
fi

nPDFsetName=$(get_opt "nPDFset" $@)
if [ -z ${nPDFset} ]; then
		nPDFsetName=EPPS21nlo_CT18Anlo_Pb208
		echo_note "nPDF set not specified. Using default: ${nPDFsetName}"
fi

cmnd_file=$(get_opt "cmnd" $@)
if [ -z ${cmnd_file} ]; then
		echo_error "Command file not specified. Exiting."
		$BASH_SOURCE --help
		exit 1
fi

output_prefix=$(get_opt "output_prefix" $@)
if [ -z ${output_prefix} ]; then
		output_prefix="$PWD/nPDFgen"
		echo_note "Output prefix not specified. Using default: ${output_prefix}"
fi

output_dir=$(get_opt "output_dir" $@)
if [ -z ${output_dir} ]; then
		_cmnd_file=$(basename ${cmnd_file})
		output_dir=${output_prefix}/${_cmnd_file}-${nPDFsetName}
		echo_note "Output directory not specified. Using default: ${output_dir}"
fi

# get the number of sets in an nPDF
nsets=$(lhapdf show ${nPDFsetName} | grep "Number of members" | awk '{print $4}')

echo_info "Number of nPDF sets in ${nPDFsetName}: ${nsets}"

echo_info "Output directory: ${output_dir}"

mkdir -pv ${output_dir}
# Create a list of commands to run in parallel
commands=()

for ((i=0; i<${nsets}; i++))
do
		cmnd_file_to_run=$(replace_nPDF ${cmnd_file} ${nPDFsetName} ${i} ${output_dir})
		_cmnd_file_to_run=$(basename ${cmnd_file_to_run})
		output_file="${output_dir}/simple_eec_${_cmnd_file_to_run%.cmnd}_output.root"
		commands+=("run_pythia ${cmnd_file_to_run} ${nev} --charged")
done

# Run the commands in parallel with a progress bar
# parallel --dry-run --bar -j 20 ::: "${commands[@]}"

nproc=$(get_opt "nproc" $@)
if [ -z ${nproc} ]; then
		nproc=5
		echo_note "Number of processes not specified. Using default: ${nproc}"
fi

parallel --bar -j ${nproc} ::: "${commands[@]}"
