#!/bin/bash

./generate_hepmc_charm.sh cmass 0.77 2>&1 | tee cmass_0.77.log &
./generate_hepmc_charm.sh cmass 1.57 2>&1 | tee cmass_1.57.log &
