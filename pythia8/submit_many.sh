#!/bin/bash

for jetpt in 40 60 80; do
	for pshower in "" "dire" "vincia"; do
		for i in {1..20}; do
			./submit_pythia8.sh /rstorage/ploskon/eec_pythia $jetpt 500 $pshower
			sleep 0.1
		done
	done
done
# End of file
