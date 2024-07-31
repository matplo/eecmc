#!/bin/bash

ncounts=1000
simple_jet_rec.py --enable-eec --ncounts ${ncounts} --nev -1 --D0mode 1 --py-hardQCDcharm --py-pthatmin 15 --jet-pt-min 15 --jet-pt-max 30 --D0-pt-min 5 --D0-pt-max 30 --py-ecm 13000 --py-monash --output pythia_eec_charm_correl.root &

ncounts=1000
simple_jet_rec.py --enable-eec --ncounts ${ncounts} --nev -1 --D0mode 1 --py-hardQCDcharm --py-pthatmin 10 --jet-pt-min 10 --jet-pt-max 15 --D0-pt-min 5 --D0-pt-max 15 --py-ecm 13000 --py-monash --output pythia_eec_charm_correl_10.root &
