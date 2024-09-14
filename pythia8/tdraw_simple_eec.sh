#!/bin/bash

. ./utils.sh

files=(
		"simple_eec_pythia_pp_5TeV_output.root"
		"simple_eec_pythia_pp_5TeV_nPDF_output.root"
		"simple_eec_pythia_pPb_5TeV_output.root"
		"simple_eec_pythia_pPb_5TeV_nPDF_output.root"
		"simple_eec_pythia_pPb_5TeV_nPDF_argantyr_output.root"
)

for fn in "${files[@]}"; do
	fn=$PWD/$fn
	if [ ! -f $fn ]; then
		echo_error "File $fn does not exist."
	else
		separator "$fn"
		../exec_tdraw_file.sh $PWD/../tdraw_eec_simple.yaml $fn
	fi
done

separator "Check the output files"

ls -ltr simple_*_h.root

separator "Done"
