#!/bin/bash

function usage()
{
	echo "Usage: $0 <nev>"
	exit 1
}

nev=$1
if [ -z $nev ]; then
	usage()
fi

./pythia8_jetreco.py --config example_config.yaml --py-pthatmin 10. --py-ecms 13000.  --jet-abs-eta-max 0.5 --jet-pt-min 10. --jet-pt-max 15. --part-abs-eta-max 0.9 --ncounts ${nev} --parton-hadron --py-hardQCDcharm --output-file charm_10_15.root &
./pythia8_jetreco.py --config example_config.yaml --py-pthatmin 15. --py-ecms 13000.  --jet-abs-eta-max 0.5 --jet-pt-min 15. --jet-pt-max 30. --part-abs-eta-max 0.9 --ncounts ${nev} --parton-hadron --py-hardQCDcharm --output-file charm_15_30.root &

./pythia8_jetreco.py --config example_config.yaml --py-pthatmin 10. --py-ecms 13000.  --jet-abs-eta-max 0.5 --jet-pt-min 10. --jet-pt-max 15. --part-abs-eta-max 0.9 --ncounts ${nev} --parton-hadron --py-hardQCD --output-file inclusive_10_15.root &
./pythia8_jetreco.py --config example_config.yaml --py-pthatmin 15. --py-ecms 13000.  --jet-abs-eta-max 0.5 --jet-pt-min 15. --jet-pt-max 30. --part-abs-eta-max 0.9 --ncounts ${nev} --parton-hadron --py-hardQCD --output-file inclusive_15_30.root &



