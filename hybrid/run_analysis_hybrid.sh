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
THIS_DIR=$(thisdir)

nev=10000

# ./analyze_hybrid.py -i HYBRID_Hadrons_Vac_5020_inc_lowpt.dat -o eec_hybrid_vac.root --nev ${nev}
# ./analyze_hybrid.py -i HYBRID_Hadrons_NoElastic_5020_inc_lowpt_05.dat --wake -o eec_hybrid_med.root --nev ${nev}
# ./analyze_hybrid.py -i HYBRID_Hadrons_NoElastic_5020_inc_lowpt_05.dat --wake --ignore-wake -o eec_hybrid_med_nw.root --nev ${nev}


#!/bin/bash

function run_analysis {
    jptmin=$1
    nev=$2
    # ./analyze_hybrid.py --log -i HYBRID_Hadrons_Vac_5020_inc_lowpt.dat -o eec_hybrid_vac_${jptmin}.root --nev ${nev} 																		--jet-min-pt ${jptmin}
    ./analyze_hybrid.py --log -i HYBRID_Hadrons_NoElastic_5020_inc_lowpt_05.dat --wake -o eec_hybrid_med_${jptmin}.root --nev ${nev} 										--jet-min-pt ${jptmin}
    ./analyze_hybrid.py --log -i HYBRID_Hadrons_NoElastic_5020_inc_lowpt_05.dat --wake --ignore-wake -o eec_hybrid_med_nw_${jptmin}.root --nev ${nev} 	--jet-min-pt ${jptmin}
		# ./analyze_hybrid.py --log --pythia --py-hardQCD --py-seed -1 --py-pthatmin ${jptmin} --nev ${nev} -o eec_pythia_${jptmin}.root 											--jet-min-pt ${jptmin} --py-ecm 5000.0
		# ./analyze_hybrid.py --log --pythia --py-hardQCD --py-seed -1 --py-pthatmin ${jptmin} --nev ${nev} -o eec_vincia_${jptmin}.root 	--py-vincia 				--jet-min-pt ${jptmin} --py-ecm 5000.0
		# ./analyze_hybrid.py --log --pythia --py-hardQCD --py-seed -1 --py-pthatmin ${jptmin} --nev ${nev} -o eec_dire_${jptmin}.root 		--py-dire 					--jet-min-pt ${jptmin} --py-ecm 5000.0
}
export -f run_analysis

parallel --eta --joblog run_hybrid.log --progress run_analysis ::: 20 40 60 80 100 ::: ${nev}

exit $?

# now drawing
# hname="hJetE2C_norm_per_bin_ch_jetR4_2040_trk10_pythia"
# hfname="${THIS_DIR}/wenqing_mc_output.root"
# rfiles=$(find . -name "*.root")
# mkdir ${THIS_DIR}/histos
# for fname in ${rfiles}
# do
# 	output="${THIS_DIR}/histos/h_$(basename ${fname})"
# 	${THIS_DIR}/draw_RL.py -i ${fname} --use-h "${hfname}:${hname}" --ncorrel 2 -o ${output}
# done