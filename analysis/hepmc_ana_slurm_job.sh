#!/bin/bash
#SBATCH --job-name=eec_sherpaD           # Job name
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks-per-node=1         # Number of tasks per node
#SBATCH --cpus-per-task=1           # Number of CPU cores per task
#SBATCH --mem=1G                    # Memory per node
#SBATCH --time=24:00:00             # Walltime
#SBATCH --partition=std  			# Job partition
#SBATCH --output={{output_dir}}/%j.out       # Output file
#SBATCH --error={{output_dir}}/%j.err        # Output file

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

cd {{output_dir}}

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

if [ ! -e {{eecmc_dir}}/eecmc.module ]; then
  echo "Error: eecmc.module not found"
  exit 1
fi
module use {{eecmc_dir}}
module load eecmc.module

# Run the application

input_file=$1
if [ -z ${input_file} ]; then
  echo "Error: input_file not found"
  exit 1
fi
# ./hepmc_D_analysis.py --config analysis_config.yaml --input ${input_file} --output ${input_file}.root --nev ${nev} 2>&1 | tee hepmc_D_analysis.log
# output_file=$(basename ${input_file} .hepmc)D.root
output_file=${input_file}_D.root
simple_jet_rec.py --enable-eec --nev ${nev} --ncounts -1 --output ${output_file} --hepmc 3 --input ${input_file} --jet-pt-min {{jet_min_pt}} --D0mode 1 --D0-pt-min 5
# output_file=$(basename ${input_file} .hepmc).root
output_file=${input_file}.root
simple_jet_rec.py --enable-eec --nev ${nev} --ncounts -1 --output ${output_file} --hepmc 3 --input ${input_file} --jet-pt-min {{jet_min_pt}}


# End of file