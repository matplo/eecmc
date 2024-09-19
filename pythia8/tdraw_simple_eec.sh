#!/bin/bash

. ./util.sh

files=(
		"simple_eec_pythia_pp_5TeV_output.root"
		"simple_eec_pythia_pp_5TeV_nPDF_output.root"
		"simple_eec_pythia_pPb_5TeV_output.root"
		"simple_eec_pythia_pPb_5TeV_nPDF_output.root"
		"simple_eec_pythia_pPb_5TeV_argantyr_output.root"
		"simple_eec_pythia_pPb_5TeV_nPDF_argantyr_output.root"
)

commands=()
for fn in "${files[@]}"; do
	fn=$PWD/$fn
	if [ ! -f $fn ]; then
		echo_error "File $fn does not exist."
	else
		separator "$fn"
		commands+=("../exec_tdraw_file.sh $PWD/../tdraw_eec_simple_smallnbin.yaml $fn ${fn%.*}_smallnbin_h.root")
		commands+=("../exec_tdraw_file.sh $PWD/../tdraw_eec_simple.yaml $fn")
	fi
done

nproc=$(get_arg "--nproc" $@)
if [ -z "${nproc}" ]; then
		nproc=6
		echo_note "Number of processes not specified. Using default: ${nproc}"
else
		echo_note "Number of processes: ${nproc}"
fi

parallel --bar -j ${nproc} ::: "${commands[@]}"

separator "Check the output files"

ls -ltr simple_*_h.root

separator "Done"
