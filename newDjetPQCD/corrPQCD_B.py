#!/usr/bin/env python

# -*- coding: utf-8 -*-

# get the hratio_pt10-15_R0.4_D0_charmNOdecays_comparison histogram from the file AnalysisResultsFinal_D0_charmNOdecays_comparison.root
# from /Users/ploskon/devel/eecmc/unfold/EEC_DKyle.root get h_10_15_logx
# divide the two histograms h_10_15_logx / hratio_pt10-15_R0.4_D0_charmNOdecays_comparison
# and save it to a file

import ROOT

# Open the ROOT files
file1 = ROOT.TFile.Open("/Users/ploskon/devel/eecmc/unfold/EEC_DKyle.root")
file2 = ROOT.TFile.Open("AnalysisResultsFinal_D0_charmNOdecays_comparison.root")

# Retrieve the histograms
hist_numerator = file1.Get("h_10_15_logx")
hist_denominator = file2.Get("hratio_pt10-15_R0.4_D0_charmNOdecays_comparison")

# Create a new histogram for the ratio
hist_ratio = hist_numerator.Clone("hist_ratio")
hist_ratio.Reset()

# Loop over the bins of the numerator histogram
for bin in range(1, hist_numerator.GetNbinsX() + 1):
    num_content = hist_numerator.GetBinContent(bin)
    x = hist_numerator.GetBinCenter(bin)
    dbin = hist_denominator.FindBin(x)
    denom_content = hist_denominator.GetBinContent(dbin)

    # Avoid division by zero
    if denom_content != 0:
        # ratio = num_content / denom_content
        ratio = num_content * denom_content
    else:
        ratio = 0

    print('x:', x, 'dbin: ', dbin, 'bin: ', bin, 'num_content: ', num_content, 'denom_content: ', denom_content, 'ratio: ', ratio)
    hist_ratio.SetBinContent(bin, ratio)

# Save the new histogram to a file
output_file = ROOT.TFile.Open("hist_ratio_output.root", "RECREATE")
hist_ratio.Write()
output_file.Close()

# Close the input files
file1.Close()
file2.Close()
