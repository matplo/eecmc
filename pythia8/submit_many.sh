#!/bin/bash

for jetpt in 20 40 60; do
    # for pshower in "" "vincia"; do
    # for pshower in "dire" "vincia"; do
    for pshower in "" "dire" "vincia"; do
		# for i in {1..100}; do
		for i in {1..100}; do
	    	./submit_pythia8.sh /rstorage/ploskon/eec_pythia $jetpt 1000 $pshower
	    	sleep 0.1
		done
    done
done
# End of file
