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
module load mysherpa

# Run the application

Sherpa -f Run.dat 2>&1 | tee sherpa_0.log
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/Process/Amegic/lib
export LD_RUN_PATH=$LD_RUN_PATH:$PWD/Process/Amegic/lib
./makelibs
nev="{{number_of_events}}"
# if nev is not an integer, then set it to 1000
# if ! [[ ${nev} =~ ^[0-9]+$ ]]; then
#   nev=1000
# fi
Sherpa -f Run.dat -e ${nev} 2>&1 | tee sherpa_1.log 
input_file=$(ls *.hepmc)
./hepmc_D_analysis.py --config analysis_config.yaml --input ${input_file} --output ${input_file}.root --nev ${nev} 2>&1 | tee hepmc_D_analysis.log

# End of file