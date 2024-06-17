#!/bin/bash

lists="sherpa_charm_hepmc_list.txt sherpa_inclusive_hepmc_list.txt"

for list in ${lists}
do
	for fn in $(cat ${list})
	do
		echo "Processing ${fn}"
		output_file=$(dirname ${fn})/$(basename ${fn} .hepmc)_ana_eec.root
		check_root_file.py ${output_file}
		if [ $? -eq 0 ]; then
			echo "File ${output_file} OK"
			continue
		else
			echo "* File ${output_file} not OK"
		fi
		echo "--> sbatch /software/users/ploskon/eecmc/analysis/hepmc_ana_slurm_job.sh ${fn}"
		sbatch /software/users/ploskon/eecmc/analysis/hepmc_ana_slurm_job.sh ${fn}
	done
done