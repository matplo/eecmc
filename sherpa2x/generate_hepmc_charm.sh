#!/bin/bash

function thisdir()
{
	SOURCE="${BASH_SOURCE[0]}"
	while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
	  DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
	  SOURCE="$(readlink "$SOURCE")"
	  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
	done
	DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
	echo ${DIR}
}
THISD=$(thisdir)

cd ${THISD}

savedir=${PWD}
nev=100000
# nev=1000

datfile_default=${THISD}/charmRun.dat
datfile=${datfile_default}
if [ "$1" == "lund" ]; then
	datfile="${THISD}/charmRunDISLundTune.dat"
	dirmod="_lund"
fi
dname_base="/rstorage/ploskon/eec_sherpa_charm_test"
mkdir -p ${dname_base}
for jetpt in 15 10
do
		dname="${dname_base}/{charm_jetpt${jetpt}${dirmod}/${nj}"
		rseed=12345
		if [ -d ${dname} ]; then
			echo "[i] directory exists - skipping - ${dname}"
		else
			mkdir ${dname}
			cd ${dname}
			yasprepl -f ${datfile} \
						--define \
						jet_min_pt="${jetpt}.0" \
						jet_eta_max=0.5 \
						jet_R=0.4 \
						p_beam_energy=2510 \
						random_seed=${rseed} \
						c_mass=${charm_mass} \
						-o Run.dat
			echo "Sherpa -f Run.dat" > ./run_sherpa.sh
			echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/Process/Amegic/lib" >> ./run_sherpa.sh
			echo "export LD_RUN_PATH=$LD_RUN_PATH:$PWD/Process/Amegic/lib" >> ./run_sherpa.sh
			echo "./makelibs" >> run_sherpa.sh
			echo "Sherpa -f Run.dat -e ${nev}" >> ./run_sherpa.sh
			#echo "„‚hepmc_jetreco_eec.py -i sherpa_LHC_jets_${jetpt}.hepmc --ncorrel 2 --hepmc 3 --output sherpa_LHC_jets_${jetpt}_eec.root"
			chmod +x ./run_sherpa.sh
			./run_sherpa.sh 2>&1 | tee sherpa.log &
			cd ${savedir}
		fi
done

cd ${savedir}
