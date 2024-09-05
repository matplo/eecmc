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

import array


def add_x_errors(graph):
    # Get the number of points in the graph
    n = graph.GetN()

    # Loop over the points in the graph
    for i in range(n):
        # Get the x and y values of the current point
        x, y = array.array("d", [0.0]), array.array("d", [0.0])
        graph.GetPoint(i, x, y)

        # Calculate the x error
        if i < n - 1:
            # Get the x value of the next point
            x_next, y_next = array.array("d", [0.0]), array.array("d", [0.0])
            graph.GetPoint(i + 1, x_next, y_next)

            # Calculate the difference between the x values
            dx = x_next[0] - x[0]

            # Set the x error to 10% of the difference
            x_error = 0.3 * dx
            # x_error = 0.5 * dx
        else:
            # For the last point, use the same x error as for the previous point
            x_error = graph.GetErrorX(i - 1)

        # Set the x error
        graph.SetPointError(
            i, x_error, x_error, graph.GetErrorYlow(i), graph.GetErrorYhigh(i)
        )

    return graph

import os
fname_10_15 = "/Users/ploskon/devel/eecmc/pQCD/fromKyle/10to15.txt"
fname_15_30 = "/Users/ploskon/devel/eecmc/pQCD/fromKyle/15to30.txt"
graph10_15 = create_graph(fname_10_15)
graph15_30 = create_graph(fname_15_30)

# Open the ROOT file
file_corr = ROOT.TFile.Open("AnalysisResultsFinal_D0_charmNOdecays_comparison.root")
# Retrieve the histograms
hist_corr_10_15 = file_corr.Get("hratio_pt10-15_R0.4_D0_charmNOdecays_comparison")
hist_corr_15_30 = file_corr.Get("hratio_pt15-30_R0.4_D0_charmNOdecays_comparison")

graph10_15_corr = create_graph(fname_10_15, hist_corr_10_15)
graph15_30_corr = create_graph(fname_15_30, hist_corr_15_30)

# Save the new histogram to a file
output_file = ROOT.TFile.Open("pQCD_Djet_corrected.root", "RECREATE")
graph10_15.SetName('pQCD_10to15')
graph10_15.SetTitle("pQCD_10to15")
graph10_15.Write()
graph10_15_corr.SetName("pQCD_10to15_corr")
graph10_15_corr.SetTitle("pQCD_10to15_corr")
graph10_15_corr.Write()

graph15_30.SetName("pQCD_15to30")
graph15_30.SetTitle("pQCD_15to30")
graph15_30.Write()
graph15_30_corr.SetName("pQCD_15to30_corr")
graph15_30_corr.SetTitle("pQCD_15to30_corr")
graph15_30_corr.Write()

output_file.Close()

# Close the input files
file_corr.Close()
