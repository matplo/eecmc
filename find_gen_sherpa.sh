#!/bin/bash

dirs="/rstorage/ploskon/eec_sherpa_charm /rstorage/ploskon/eec_sherpa_inclusive"

regenerate=$1
if [ "x$regenerate" == "xregenerate" ]; then
	echo "Regenerating lists..."
	rm sherpa_charm_hepmc_list.txt sherpa_inclusive_hepmc_list.txt
	touch sherpa_charm_hepmc_list.txt sherpa_inclusive_hepmc_list.txt
	for dn in $dirs 
	do
		if [[ $dn == *"sherpa_charm"* ]]; then
			fname="sherpa_charm_hepmc_list.txt"
			fname_root="sherpa_charm_root_list.txt"
		fi
		if [[ $dn == *"sherpa_inclusive"* ]]; then
			fname="sherpa_inclusive_hepmc_list.txt"
			fname_root="sherpa_inclusive_root_list.txt"
		fi
		echo "Processing ${dn} - hepmc..."
		find $dn -name "*.hepmc" | tee -a ${fname}
		echo "Processing ${dn} - root..."
		find $dn -name "*_ana_eec.root"  | tee -a | tee -a ${fname_root}
	done
fi

list=$(cat sherpa_charm_hepmc_list.txt sherpa_inclusive_hepmc_list.txt)

for input_file in $list
do
	output_file="list_"
	# setup jet pt cuts
	if [[ ${input_file} == *"5020"* ]]; then
		ecm="5TeV"
	fi

	if [[ ${input_file} == *"13000"* ]]; then
		ecm="13TeV"
	fi

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

	fragm=""
	if [[ ${input_file} == *"lund"* ]]; then
		fragm="lund"
	fi

	flav=""
	# check if charm or inclusive
	if [[ ${input_file} == *"charm"* ]]; then
		flav="charm"
	else
		flav="inclusive"
	fi
	jet_setting="jetpt_${jet_min_pt}"
	if [ ! -z $fragm ]; then
		jet_setting="${fragm}_jetpt_${jet_min_pt}"
	fi
	output_file="${output_file}_${flav}_${ecm}_${jet_setting}.txt"
	echo ${input_file} | tee -a ${output_file}
done
