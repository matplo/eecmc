# ../analysis/simple_jet_rec.py --enable-eec --ncounts 100000 --nev 10000 --D0mode 1 --debug

# ../analysis/simple_jet_rec.py --enable-eec --ncounts 100000 --nev 10000 --D0mode 1 --debug --py-hardQCDcharm

# ../analysis/simple_jet_rec.py --enable-eec --ncounts 1000 --nev 1000 --D0mode 1 --py-hardQCDcharm $@

# ../analysis/simple_jet_rec.py --enable-eec --ncounts 1000 --nev 1000 --D0mode 1 --py-hardQCDcharm --py-pthatmin 15 --jet-pt-min 15 --jet-pt-max 30 --D0-pt-min 5 --D0-pt-max 30 $@

../analysis/simple_jet_rec.py --enable-eec --ncounts 1000 --nev -1 --D0mode 1 --py-hardQCDcharm --py-pthatmin 15 --jet-pt-min 15 --jet-pt-max 30 --D0-pt-min 5 --D0-pt-max 30 --nev 100000 --py-charm-mass 0 --output pythia_charm0.root
