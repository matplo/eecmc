#!/bin/bash

# ncounts=30000
ncounts=20000

./eec_correl.py --enable-eec --ncounts ${ncounts} --nev -1 --py-hardQCD --py-pthatmin 10 --jet-pt-min 10 --jet-pt-max 15 --py-ecm 13000 --py-monash --output pythia_eec_inclusive_correl_10.root &

./eec_correl.py --enable-eec --ncounts ${ncounts} --nev -1 --py-hardQCD --py-pthatmin 15 --jet-pt-min 15 --jet-pt-max 30 --py-ecm 13000 --py-monash --output pythia_eec_inclusive_correl_15.root &

# D0 modes:
# 0 - no D0 selection
# 1 - D0 selection
# 2 - D0 selection and remove D0 daughters
# 3 - D0 selection Kpi decay
# 4 - D0 selection Kpi decay and remove D0 daughters

ncounts=5000

for mode in 0 1 2 3 4
do
	./eec_correl.py --enable-eec --ncounts ${ncounts} --nev -1 --D0mode ${mode} --py-hardQCDcharm --py-pthatmin 10 --jet-pt-min 10 --jet-pt-max 15 --D0-pt-min 5 --D0-pt-max 15 --py-ecm 13000 --py-monash --output pythia_eec_charm_correl_10_D${mode}.root &
	./eec_correl.py --enable-eec --ncounts ${ncounts} --nev -1 --D0mode ${mode} --py-hardQCDcharm --py-pthatmin 15 --jet-pt-min 15 --jet-pt-max 30 --D0-pt-min 5 --D0-pt-max 30 --py-ecm 13000 --py-monash --output pythia_eec_charm_correl_10_D${mode}.root &
done

#D0 mode 2 is selecting events with D0 AND replacing the two daughters with D0
#ncounts=2000
#./eec_correl.py --enable-eec --ncounts ${ncounts} --nev -1 --D0mode 2 --py-hardQCDcharm --py-pthatmin 10 --jet-pt-min 10 --jet-pt-max 15 --D0-pt-min 5 --D0-pt-max 15 --py-ecm 13000 --py-monash --output pythia_eec_charm_correl_10_D2_debug.root &
