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
nev=1000

if [ "$1" == "clean" ]; then
	rm -rf pythia8_jetpt*
	exit 0
fi


# for jetpt in 20
for jetpt in 40 60 80
do
	dname="pythia8_jetpt${jetpt}${dirmod}"
	if [ -d ${dname} ]; then
		echo "[i] directory exists - skipping - ${dname}"
	else
		mkdir ${dname}
		cd ${dname}
		echo "$THISD/pythia8_jetreco_eec.py --py-pthatmin $jetpt --nev $nev" >> ./run_pythia8.sh
		echo "$THISD/pythia8_jetreco_eec.py --py-pthatmin $jetpt --nev $nev --py-vincia" >> ./run_pythia8.sh
		echo "$THISD/pythia8_jetreco_eec.py --py-pthatmin $jetpt --nev $nev --py-dire" >> ./run_pythia8.sh
		chmod +x ./run_pythia8.sh
		./run_pythia8.sh 2>&1 | tee pythia8.log &
		cd ${savedir}
	fi
done

cd ${savedir}
