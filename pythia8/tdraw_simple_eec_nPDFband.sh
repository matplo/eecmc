#!/bin/bash

. ./util.sh

input_directory=$1
if [ -z ${input_directory} ]; then
	echo_error "Input directory not provided."
	echo_info "Usage: $0 <input_directory> [--nproc=<nproc>]"
	exit 1
fi

abspath_input_directory=$(realpath ${input_directory})
if [ ! -d ${abspath_input_directory} ]; then
	echo_error "Input directory ${abspath_input_directory} does not exist."
	exit 1
fi

separator "Input directory: ${abspath_input_directory}"
input_directory=${abspath_input_directory}
files=$(find ${input_directory} -name "*.root" | grep -v ".*_h.root")

echo ${files}

commands=()
for fn in $files; do
	fn=$fn
	if [ ! -f $fn ]; then
		echo_error "File $fn does not exist."
	else
		separator "$fn"
		# ../exec_tdraw_file.sh $PWD/../tdraw_eec_simple.yaml $fn
		commands+=("../exec_tdraw_file.sh $PWD/../tdraw_eec_simple.yaml $fn")
	fi
done

nproc=$(get_arg "--nproc" $@)
if [ -z "${nproc}" ]; then
		nproc=5
		echo_note "Number of processes not specified. Using default: ${nproc}"
else
		echo_note "Number of processes: ${nproc}"
fi

parallel --bar -j ${nproc} ::: "${commands[@]}"

separator "Check the output files"

ls -ltr ${input_directory}/*_h.root | wc -l

separator "Done"
