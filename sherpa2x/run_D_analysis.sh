#!/bin/bash

# ./hepmc_D_analysis.py --config analysis_config_D_jetpt10.yaml &
# ./hepmc_D_analysis.py --config analysis_config_D_jetpt15.yaml &

# ./hepmc_D_analysis.py --config analysis_config_D_jetpt10_lund.yaml &
# ./hepmc_D_analysis.py --config analysis_config_D_jetpt15_lund.yaml &

./hepmc_D_analysis.py --config analysis_config_D_jetpt15.yaml --output nc100_std_mass.root --ncounts 100 --nev 1000000
./hepmc_D_analysis.py --config analysis_config_D_jetpt15.yaml --output nc100_high_mass.root --input charm_jetpt15_cmass_1.57/sherpa_LHC_jets_15.0.hepmc --ncounts 100 --nev 1000000
./hepmc_D_analysis.py --config analysis_config_D_jetpt15.yaml --output nc100_low_mass.root --input charm_jetpt15_cmass_0.77/sherpa_LHC_jets_15.0.hepmc --ncounts 100 --nev 1000000
