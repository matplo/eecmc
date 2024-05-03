# ./hepmc_D_test.py -i sherpa_LHC_jets_30.0.hepmc --hepmc 3 --nev 10000 --jet-min-pt 7 --jet-max-pt 10 --D0-min-pt 3 --D0-max-pt 10
./hepmc_D_test.py -i sherpa_LHC_jets_30.0.hepmc --hepmc 3 --nev 10000 --jet-min-pt 30 --jet-max-pt 100 --D0-min-pt 0 --D0-max-pt 100 --nev -1 $@
