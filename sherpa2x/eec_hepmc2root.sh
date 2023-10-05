#!/bin/bash

savedir=${PWD}
nev=100000

if [ "$1" == "lund" ]; then
	dirmod="_lund"
fi

for jetpt in 20 40 60
do
	dname="jetpt${jetpt}${dirmod}"
	if [ -d ${dname} ]; then
		#cd ${dname}
		input_file=${dname}/sherpa_LHC_jets_${jetpt}.hepmc
		if [ "x${dirmod}" == "x_lund" ]; then
			input_file=$(find ${dname} -name "sherpa_LHC_jets_*.hepmc")
		fi
		#${HEPPYY_PYTHON_LIB}/heppyy/sandbox/hepmc_jetreco_eec.py -i ${input_file} --hepmc 3 --nev ${nev} -o "eec_sherpa_LHC_jetpt${jetpt}.root"
		./hepmc_jetreco_eec.py -i ${input_file} --hepmc 3 -o "${dname}/eec_sherpa_LHC_jetpt${jetpt}.root" &
		#cd ${savedir}
	fi
done

