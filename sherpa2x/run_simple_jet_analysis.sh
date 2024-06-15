#!/bin/bash

input=$1
if [ -z ${input} ]; then
	echo "Error: input file not provided"
	exit 1
fi

output=$1.root

# this assumes the eecmc.module is loaded
simple_jet_rec.py --enable-eec --ncounts -1 --D0mode 1 --output ${output} --hepmc 3 --input ${input} --jet-pt-min 15 --D0-pt-min 5
