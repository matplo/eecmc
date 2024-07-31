#!/bin/bash

files=$(ls hadded_list__*_ana_eec.txt.root)
for fn in ${files}
do
		echo ${fn}
		ptmin=15
		ptmax=30
		if [[ ${fn} == *"5TeV"* ]]; then
				ptmin=10
				ptmax=20
		elif [[ ${fn} == *"13TeV"* ]]; then
				ptmin=15
				ptmax=30
				if [[ ${fn} == *"jetpt_10"* ]]; then
						ptmin=10
						ptmax=15
				fi
		fi
		# tdraw_eec.yaml
		# tdraw_eec_pt_template.yaml
		yaml_tdraw_file="tdraw_eec_tmp.yaml"

		yaml_tdraw_template="tdraw_eec_template.yaml"
		lre ${yaml_tdraw_template} --define jet_min_pt=${ptmin} jet_max_pt=${ptmax} > ${yaml_tdraw_file}
		fn_out=${fn%.*}_h.root
		./exec_tdraw_file.sh ${yaml_tdraw_file} ${fn} ${fn_out}
		if [ $? -ne 0 ]; then
				echo "Error in ${yaml_tdraw_file}"
				exit 1
		fi

		yaml_tdraw_template="tdraw_eec_hbin_template.yaml"
		lre ${yaml_tdraw_template} --define jet_min_pt=${ptmin} jet_max_pt=${ptmax} > ${yaml_tdraw_file}
		fn_out=${fn%.*}_hbin.root
		./exec_tdraw_file.sh ${yaml_tdraw_file} ${fn} ${fn_out}
		if [ $? -ne 0 ]; then
				echo "Error in ${yaml_tdraw_file}"
				exit 1
		fi

done
