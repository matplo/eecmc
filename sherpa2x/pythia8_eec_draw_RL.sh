#!/bin/bash

savedir=${PWD}

if [ "$1" == "lund" ]; then
	dirmod="_lund"
fi

declare -a cmnds=()
for jetpt in 40 60 80
do
	dname="pythia8_jetpt${jetpt}"
	if [ -d ${dname} ]; then
		cd ${dname}
		echo "${PWD}"
		for dirmod in "" "_vincia" "_dire"
		do 
			fin="eec_pythia8${dirmod}.root"
			echo "input file: ${fin}"
			output="h_eec_pythia8${dirmod}.root"
			# cmnd="${savedir}/draw_RL.py -i ${PWD}/${fin} --nbins 80 --ncorrel 3 -o ${PWD}/${output}" 
			[ "${jetpt}" == "20" ] && hname="hJetE2C_norm_per_bin_ch_jetR4_2040_trk10_pythia"
			[ "${jetpt}" == "40" ] && hname="hJetE2C_norm_per_bin_ch_jetR4_4060_trk10_pythia"
			[ "${jetpt}" == "60" ] && hname="hJetE2C_norm_per_bin_ch_jetR4_6080_trk10_pythia"
			hfname="${savedir}/wenqing_mc_output.root"
			if [ -f "${fin}" ]; then
				cmnd="${savedir}/draw_RL.py -i ${PWD}/${fin} --use-h "${hfname}:${hname}" --ncorrel 2 -o ${PWD}/${output}" 
				echo ${cmnd} | tee -a ${savedir}/cmnds.txt
				cmnds[${#cmnds[@]}]=${cmnd}
			fi
		done
		cd ${savedir}
	fi
done

length=${#cmnds[@]}
for (( j=0; j<${length}; j++ ));
do
  printf "Command %d is %s\n" $j "${cmnds[$j]}"
  #parallel {1} ::: "${cmnds[$j]}"
  ${cmnds[$j]} &
done

