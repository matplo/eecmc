#!/usr/bin/env python

# -*- coding: utf-8 -*-

# get the hratio_pt10-15_R0.4_D0_charmNOdecays_comparison histogram from the file AnalysisResultsFinal_D0_charmNOdecays_comparison.root
# from /Users/ploskon/devel/eecmc/unfold/EEC_DKyle.root get h_10_15_logx
# divide the two histograms h_10_15_logx / hratio_pt10-15_R0.4_D0_charmNOdecays_comparison
# and save it to a file

import ROOT


def create_graph(filename, hcorrection=None):
    # Open the file
    with open(filename, "r") as f:
        lines = f.readlines()

    # Create a TGraphAsymmErrors object
    graph = ROOT.TGraphAsymmErrors(len(lines))

    RLs = []
    # Loop over the lines in the file
    for i, line in enumerate(lines):
        # Split the line into columns
        columns = line.split()

        # Get the values from the columns
        RL = float(columns[0])
        RLs.append(RL)
        min_theory = float(columns[1])
        central_theory = float(columns[2])
        max_theory = float(columns[3])

        if hcorrection is not None:
            xbin = hcorrection.FindBin(RL)
            min_theory = min_theory * hcorrection.GetBinContent(xbin)
            central_theory = central_theory * hcorrection.GetBinContent(xbin)
            max_theory = max_theory * hcorrection.GetBinContent(xbin)

        # Calculate the y errors
        y_error_low = central_theory - min_theory
        y_error_high = max_theory - central_theory

        # Set the point and its errors
        graph.SetPoint(i, RL, central_theory)
        graph.SetPointError(i, 0, 0, y_error_low, y_error_high)  # no error in x

    return graph

import os
fname_10_15 = "/Users/ploskon/devel/eecmc/pQCD/fromKyle/10to15.txt"
fname_15_30 = "/Users/ploskon/devel/eecmc/pQCD/fromKyle/15to30.txt"
graph10_15 = create_graph(fname_10_15)
graph15_30 = create_graph(fname_15_30)

# Open the ROOT files
file1 = ROOT.TFile.Open("/Users/ploskon/devel/eecmc/unfold/EEC_DKyle.root")
file2 = ROOT.TFile.Open("AnalysisResultsFinal_D0_charmNOdecays_comparison.root")

# Retrieve the histograms
hist_numerator = file1.Get("h_10_15_logx")
hist_denominator = file2.Get("hratio_pt10-15_R0.4_D0_charmNOdecays_comparison")

graph10_15_corr = create_graph(fname_10_15, hist_denominator)
graph15_30_corr = create_graph(fname_15_30, hist_denominator)

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
graph10_15.SetName('pQCD_10to15')
graph10_15.Write()
graph10_15_corr.SetName("pQCD_10to15_corr")
graph10_15_corr.Write()
output_file.Close()

# Close the input files
file1.Close()
file2.Close()
