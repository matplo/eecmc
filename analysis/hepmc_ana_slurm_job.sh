#!/bin/bash
#SBATCH --job-name=eec__ana           # Job name
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks-per-node=1         # Number of tasks per node
#SBATCH --cpus-per-task=1           # Number of CPU cores per task
#SBATCH --mem=1G                    # Memory per node
#SBATCH --time=04:00:00             # Walltime
#SBATCH --partition=quick  			# Job partition
#SBATCH --output=/tmp/%j.out       # Output file
#SBATCH --error=/tmp/%j.err        # Output file

input_file=$1
if [ ! -f ${input_file} ]; then
  echo "Error: input_file not found"
  exit 1
fi

EECMC_DIR=/software/users/ploskon/eecmc
if [ -z ${input_file} ]; then
  echo "Error: input_file not found"
  return 1
fi
echo "Running on ${input_file}"
# setup jet pt cuts
if [[ ${input_file} == *"jetpt10"* ]]; then
  jet_min_pt=10
  jet_max_pt=15
  if [[ ${input_file} == *"5020"* ]]; then
    jet_max_pt=20
  fi
else
  jet_min_pt=15
  jet_max_pt=30
fi
charm_setting=""
# check if charm or inclusive
if [[ ${input_file} == *"charm"* ]]; then
  echo "Running on charm"
  charm_setting="--D0mode 1 --D0-pt-min 5 --D0-pt-max ${jet_max_pt}"
else
  echo "Running on inclusive"
fi
jet_setting="--jet-pt-min ${jet_min_pt} --jet-pt-max ${jet_max_pt}"
echo "Jet setting: ${jet_setting}"
echo "Charm setting: ${charm_setting}"

# setup the environment
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
module load bundle/hepbase
module load heppyy/current

if [ ! -e $EECMC_DIR/eecmc.module ]; then
  echo "Error: $EECMC_DIR/eecmc.module not found"
  exit 1
fi
module use $EECMC_DIR
module load eecmc.module

output_file=$(dirname ${input_file})/$(basename ${input_file} .hepmc)_ana_eec.root
script_file=$(dirname ${input_file})/$(basename ${input_file} .hepmc)_ana_eec.sh
log_file=$(dirname ${input_file})/$(basename ${input_file} .hepmc)_ana_eec.log
nev=10000
echo "simple_jet_rec.py --enable-eec --nev ${nev} --output ${output_file} --hepmc 3 --input ${input_file} ${jet_setting} ${charm_setting}" | tee ${script_file}
chmod +x ${script_file}
# ${script_file} 2>&1 | tee ${log_file}
simple_jet_rec.py --enable-eec --nev ${nev} --output ${output_file} --hepmc 3 --input ${input_file} ${jet_setting} ${charm_setting} 2>&1 | tee ${log_file}


# End of file