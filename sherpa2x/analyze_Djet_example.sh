# ./hepmc_D_test.py -i sherpa_LHC_jets_30.0.hepmc --hepmc 3 --nev 10000 --jet-min-pt 7 --jet-max-pt 10 --D0-min-pt 3 --D0-max-pt 10
finput=$1
if [ -z $finput ]; then
	echo "Usage: $0 <input.hepmc>"
	exit 1
fi
shift 1
# ./hepmc_D_test.py -i $finput --hepmc 3 --nev 10000 --jet-min-pt 30 --jet-max-pt 100 --D0-min-pt 0 --D0-max-pt 100 --nev -1 $@
./hepmc_D_test.py -i $finput --hepmc 3 --nev 10000 --jet-min-pt 15 --jet-max-pt 30 --D0-min-pt 5 --D0-max-pt 6 --nev -1 -o 15-30-D-5-6.root $@
