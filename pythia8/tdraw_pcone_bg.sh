#!/bin/bash

# files="pcone_bg_pp.root pcone_bg_pPb_argantyr.root"
files="pcone_bg_pPb.root"

for fn in $files; do
	fn=$PWD/$fn
	if [ ! -f $fn ]; then
		echo_error "File $fn does not exist."
	else
		separator "$fn"
		../exec_tdraw_file.sh $PWD/../tdraw_eec_simple_smallnbin.yaml $fn ${fn%.*}_smallnbin_h.root
		../exec_tdraw_file.sh $PWD/../tdraw_eec_simple.yaml $fn
	fi
done
