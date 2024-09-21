#!/bin/bash

. ./util.sh

nev=100000
commands=()
commands+=("./pythia8_simple_eec.py --py-cmnd pythia_pp_5TeV.cmnd --jet-pt-min 20. --charged --py-pthatmin 20. --py-seed 123456 -o pcone_bg_pp.root --nev ${nev} 2>&1 | tee pcone_bg_pp.log")
commands+=("./pythia8_simple_eec.py --py-cmnd pythia_pPb_5TeV_nPDF_argantyr.cmnd --jet-pt-min 20. --charged --py-pthatmin 20. --py-seed 123456 -o pcone_bg_pPb_argantyr.root --nev ${nev} 2>&1 | tee pcone_bg_pPb_argantyr.log")

nproc=$(get_arg "--nproc" $@)
if [ -z "${nproc}" ]; then
		nproc=2
		echo_note "Number of processes not specified. Using default: ${nproc}"
else
		echo_note "Number of processes: ${nproc}"
fi

parallel --bar -j ${nproc} ::: "${commands[@]}"
separator "Done"

