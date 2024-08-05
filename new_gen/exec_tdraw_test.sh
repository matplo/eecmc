#!/bin/bash

#  KEY: TNtuple  tn_jet_parton;1 tn_jet
#  KEY: TNtuple  tn_jet_parton_eec;1     tn_eec
#  KEY: TNtuple  tn_jet_ch_parton_eec;1  tn_eec_ch
#  KEY: TNtuple  tn_jet_hadron;1 tn_jet
#  KEY: TNtuple  tn_jet_hadron_eec;1     tn_eec
#  KEY: TNtuple  tn_jet_ch_hadron_eec;1  tn_eec_ch
#  KEY: TNtuple  tn_jet_charged;1        tn_jet
#  KEY: TNtuple  tn_jet_charged_eec;1    tn_eec
#  KEY: TNtuple  tn_jet_ch_charged_eec;1 tn_eec_ch
#  KEY: TNtuple  tn_jet_D0;1     tn_jet
#  KEY: TNtuple  tn_jet_D0_eec;1 tn_eec
#  KEY: TNtuple  tn_jet_ch_D0_eec;1      tn_eec_ch
#  KEY: TNtuple  tn_jet_D0Kpi;1  tn_jet
#  KEY: TNtuple  tn_jet_D0Kpi_eec;1      tn_eec
#  KEY: TNtuple  tn_jet_ch_D0Kpi_eec;1   tn_eec_ch

file_input=$1
if [ -z "$file_input" ]; then
		echo "[e] file_input missing"
		exit 1
fi

file_output=$2
if [ -z "$file_output" ]; then
		file_output=${file_input%.*}_h.root
fi

for tn_jet in tn_jet_parton tn_jet_hadron tn_jet_charged tn_jet_D0 tn_jet_D0Kpi
do
	foutput=${file_output%.*}_${tn_jet}.root
	yaml_file=${foutput%.*}_tdraw_conf.yaml
	lre tdraw_conf.yaml --define input="$file_input" output="$foutput" jet_min_pt=10 jet_max_pt=15 tn_jet=${tn_jet} > ${yaml_file}
	../draw_from_yaml.py -c ${yaml_file}
done