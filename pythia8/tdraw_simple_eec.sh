#!/bin/bash

../exec_tdraw_file.sh $PWD/../tdraw_eec_simple.yaml $PWD/simple_eec_pythia_pp_5TeV_output.root 
../exec_tdraw_file.sh $PWD/../tdraw_eec_simple.yaml $PWD/simple_eec_pythia_pp_5TeV_nPDF_output.root 
../exec_tdraw_file.sh $PWD/../tdraw_eec_simple.yaml $PWD/simple_eec_pythia_pPb_5TeV_output.root 
../exec_tdraw_file.sh $PWD/../tdraw_eec_simple.yaml $PWD/simple_eec_pythia_pPb_5TeV_nPDF_output.root 

ls -ltr simple_*_h.root