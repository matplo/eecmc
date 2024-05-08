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
module avail
module load yasp
module load mysherpa

# ./hepmc_D_analysis.py --config analysis_config_D_jetpt10.yaml &
./hepmc_D_analysis.py --config analysis_config_D_jetpt15.yaml &

# ./hepmc_D_analysis.py --config analysis_config_D_jetpt10_lund.yaml &
./hepmc_D_analysis.py --config analysis_config_D_jetpt15_lund.yaml &

# ./hepmc_D_analysis.py --config analysis_config_D_jetpt15.yaml --output nc100_std_mass.root --ncounts 100 --nev 1000000
# ./hepmc_D_analysis.py --config analysis_config_D_jetpt15.yaml --output nc100_high_mass.root --input charm_jetpt15_cmass_1.57/sherpa_LHC_jets_15.0.hepmc --ncounts 100 --nev 1000000
# ./hepmc_D_analysis.py --config analysis_config_D_jetpt15.yaml --output nc100_low_mass.root --input charm_jetpt15_cmass_0.77/sherpa_LHC_jets_15.0.hepmc --ncounts 100 --nev 1000000
