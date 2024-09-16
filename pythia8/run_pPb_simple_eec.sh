#!/bin/bash

run_pythia() {
    local cmnd_file=$1
    local nev=$2
    local charged=$3
    local output_file="simple_eec_${cmnd_file%.cmnd}_output.root"
    
    if [ ! -f $cmnd_file ]; then
        echo "File $cmnd_file does not exist."
        exit 1
    fi
    
    ./pythia8_simple_eec.py --py-cmnd ${cmnd_file} --nev ${nev} --jet-pt-min 20. ${charged} --py-pthatmin 20. -o ${output_file} 2>&1 | tee ${output_file%.root}.log
}

export -f run_pythia

nev=100000
charged="--charged"
# charged=""

#    "pythia_pPb_8TeV.cmnd"

# Define the command files
cmnd_files=(
    "pythia_pPb_5TeV.cmnd"
    "pythia_pPb_5TeV_nPDF.cmnd"
    "pythia_pPb_5TeV_argantyr.cmnd"
	"pythia_pPb_5TeV_nPDF_argantyr.cmnd"
)

# Create a list of commands to run in parallel
commands=()
for cmnd_file in "${cmnd_files[@]}"; do
    commands+=("run_pythia ${cmnd_file} ${nev} ${charged}")
    # Replace _pPb_ with _pp_
    pp_cmnd_file=${cmnd_file/pPb/pp}
    commands+=("run_pythia ${pp_cmnd_file} ${nev} ${charged}")
done

# Run the commands in parallel with a progress bar
# parallel --dry-run --bar ::: "${commands[@]}"
parallel --bar ::: "${commands[@]}"