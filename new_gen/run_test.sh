#!/bin/bash

#--debug

ncounts=100000

./eec_correl.py --py-monash --py-pthatmin 10 --py-hardQCDcharm --ncounts $ncounts --nev -1 --jet-pt-min 10 --jet-pt-max 15 --D0-pt-min 5 --D0-pt-max 15 --jet-abs-eta-max 0.5 -o test_10_100k.root $@ &
./eec_correl.py --py-monash --py-pthatmin 15 --py-hardQCDcharm --ncounts $ncounts --nev -1 --jet-pt-min 15 --jet-pt-max 30 --D0-pt-min 5 --D0-pt-max 30 --jet-abs-eta-max 0.5 -o test_15_100k.root $@ &
