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
#module avail
module load yasp
module load bundle/hepbase
# module load bundle/sherpa2x
module load heppyy/current

lre module.template --define EECMCDIR=$THISD > eecmc.module

cd ${savedir}
# End of make_module.sh
