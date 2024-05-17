#!/usr/bin/env python

from uproot_tree import UprootTree
t = UprootTree()

t.add("inclusive_10_15_jfch.root", "tn_events_jfch", 'ev')
t.add_cut("ev", "ev_njets", 1, 1000)

t.add('inclusive_10_15_jfch_jet_full.root', 'tn_jet_jet_full', 'jet')

# t.add('inclusive_10_15_jfch_jet_full.root', 'tn_parts_jet_full', 'part')

# t.add("inclusive_10_15_jfch.root", "tn_correl_jfch", 'corr')
# t.add_cut("corr_pt", 0, 1000)

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

# print(t.get_df())

for group_name, group_df in t.get_df():
    df = group_df
    # print_group(group_name, group_df)
# how to fill a histogram:
import ROOT
import numpy as np
import array

# Assuming df is your DataFrame and 'column' is the name of the column
data = []
try:
    while True:
        dname, df = t.next()
        if df is None:
            break
        # print(df)
        _ = [data.append(x) for x in df['jet_pt'].values]
except StopIteration:
    pass

data = np.concatenate(t.get_unique(['jet_pt']).values).tolist()
print(type(data))
# print(data)
# Create a ROOT histogram
hist = ROOT.TH1F('hist', 'Histogram from pandas DataFrame', 100, min(data), max(data) + 1)
# hist = ROOT.TH1F('hist', 'Histogram from pandas DataFrame', 100, data.min(), data.max() + 1)

# Convert numpy arrays to C++ compatible arrays
data_array = array.array('d', data)
weights_array = array.array('d', np.ones(len(data)))

# Fill the histogram
hist.FillN(len(data), data_array, weights_array)

fx = ROOT.TFile('inclusive_10_15_jfch_jet_full.root')
tx = fx.Get('tn_jet_jet_full')
hist2 = ROOT.TH1F('hist2', 'Histogram from pandas DataFrame', 100, min(data), max(data) + 1)
tx.Draw('pt>>hist2')

print("compare n entries:", hist.GetEntries())
print("compare n entries:", hist2.GetEntries())

# Draw the histogram
canvas = ROOT.TCanvas('c', 'c', 800, 600)
hist2.SetLineColor(ROOT.kYellow)
hist2.SetLineWidth(10)
hist2.Draw()
hist.Draw('same')
canvas.Draw()
canvas.SaveAs('hist.png')