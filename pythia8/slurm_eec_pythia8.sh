#!/bin/bash
#SBATCH --job-name=eec_pythia8           # Job name
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks-per-node=1         # Number of tasks per node
#SBATCH --cpus-per-task=1           # Number of CPU cores per task
#SBATCH --mem=1G                    # Memory per node
#SBATCH --time=02:00:00             # Walltime
#SBATCH --partition=quick  			# Job partition
#SBATCH --output=slurm_%j.out       # Output file
#SBATCH --error=slurm_%j.err        # Output file

# Load any necessary modules
# SYS_YASP_DIR=$(yaspenv.sh yasp -q feature yasp_dir 2>&1 | tail -n 1)
SYS_YASP_DIR={{SYS_YASP_DIR}}
source ${SYS_YASP_DIR}/venvyasp/bin/activate
module use ${SYS_YASP_DIR}/software/modules
module load yasp myheppyy

# Run the application

cd {{OUTPUT_DIR}}
$HEPPYY_DEV/heppyy/example/test_yaspcppyy_pythia_fastjet_simple_load.py --nev 10

# End of file