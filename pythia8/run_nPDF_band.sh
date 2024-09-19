#!/bin/bash

savedir=$PWD

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

yaspenv_shell=$(which yaspenv.sh)
if [ -z ${yaspenv_shell} ]; then
  echo "Error: yaspenv.sh not found"
  exit 1
fi

SYS_YASP_DIR=$(${yaspenv_shell} yasp -q feature yasp_dir 2>&1 | tail -n 1)
if [ -z ${SYS_YASP_DIR} ]; then
  echo "Error: SYS_YASP_DIR not found"
  exit 1
fi

echo "[i] using ${SYS_YASP_DIR}"
source ${SYS_YASP_DIR}/venvyasp/bin/activate
module use ${SYS_YASP_DIR}/software/modules
module load yasp
module load bundle/hepbase
module load heppyy/current

EEC_DIR_TMP=$(dirname ${THISD})
if [ ! -e ${EEC_DIR_TMP}/eecmc.module ]; then
  echo "Error: ${EEC_DIR_TMP}/eecmc.module not found"
  exit 1
fi
module use ${EEC_DIR_TMP}
module load eecmc.module
module list

# this was just a setup...
nev=100000
./run_pPb_nPDF_uncert.sh --cmnd=pythia_pPb_5TeV_nPDF.cmnd --nev=${nev} --nproc=48
./run_pPb_nPDF_uncert.sh --cmnd=pythia_pPb_5TeV_nPDF_argantyr.cmnd --nev=${nev} --nproc=48

cd ${savedir}