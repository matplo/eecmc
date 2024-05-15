#!/usr/bin/env python

from uproot_tree import UprootTree
t = UprootTree()
t.add("test_jfch.root", "tn_correl_jfch", 'corr')
t.add_cut("corr_pt", 20, 1000) # ambigous when more than one pt column

t.add("test_jfch.root", "tn_events_jfch", 'ev')
t.add_cut("ev_njets", 2, 1000)

t.add('test_jfch_jet_full.root', 'tn_jet_jet_full', 'jet')

# t.add('test_jfch_jet_full.root', 'tn_parts_jet_full', 'part')

t.add_group_by('nev', 'xsec', 'ev_weight')

import ROOT
h = ROOT.TH1F("pt", "pt", 1000, 0, 1000)

def print_group(group_name, group_df, separator=False):
    print('group name:', group_name)
    print('group data:\n', group_df.to_string(index=True))
    if separator:
        print('--')
            
# try:
# 	nev = 0
# 	while True:
# 		# x = t.next()
# 		group_name, group_df = t.next()
# 		print_group(group_name, group_df, True)
# 		nev += 1
# except StopIteration:
# 	pass

# t.get_df().to_root("test_jfch.root", "tn_correl_jfch", mode="update")

print(t.get_df())

for group_name, group_df in t.get_df():
    df = group_df
    print_group(group_name, group_df)

# how to fill a histogram:
import ROOT
import numpy as np
import array

# Assuming df is your DataFrame and 'column' is the name of the column
data = df['jet_pt'].values

# Create a ROOT histogram
hist = ROOT.TH1F('hist', 'Histogram from pandas DataFrame', 100, data.min(), data.max()+1)

# Convert numpy arrays to C++ compatible arrays
data_array = array.array('d', data)
weights_array = array.array('d', np.ones(len(data)))

# Fill the histogram
hist.FillN(len(data), data_array, weights_array)

# Draw the histogram
canvas = ROOT.TCanvas('c', 'c', 800, 600)
hist.Draw()
canvas.Draw()
canvas.SaveAs('hist.png')