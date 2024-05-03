#!/bin/bash
#SBATCH --job-name=eec_pythia8           # Job name
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks-per-node=1         # Number of tasks per node
#SBATCH --cpus-per-task=1           # Number of CPU cores per task
#SBATCH --mem=1G                    # Memory per node
#SBATCH --time=24:00:00             # Walltime
#SBATCH --partition=std  			# Job partition
#SBATCH --output={{OUTPUT_DIR}}/slurm_%j.out       # Output file
#SBATCH --error={{OUTPUT_DIR}}/slurm_%j.err        # Output file

# Load any necessary modules
# SYS_YASP_DIR=$(yaspenv.sh yasp -q feature yasp_dir 2>&1 | tail -n 1)
source {{SYS_YASP_DIR}}/venvyasp/bin/activate
module use {{SYS_YASP_DIR}}/software/modules
module load yasp myheppyy
module avail
module list
# Run the application

cd {{OUTPUT_DIR}}
Sherpa -f Run.dat
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/Process/Amegic/lib
export LD_RUN_PATH=$LD_RUN_PATH:$PWD/Process/Amegic/lib
./makelibs
Sherpa -f Run.dat -e ${nev}
chmod +x ./run_sherpa.sh
./run_sherpa.sh 2>&1 | tee sherpa.log &
input_file=$(ls *.hepmc)
./hepmc_D_test.py -i ${input_file} --hepmc 3 --nev -1 --config analysis_config.yaml

# End of file