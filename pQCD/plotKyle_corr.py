#!/usr/bin/env python

import ROOT

def create_graph(filename):
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
        x, y = array.array('d', [0.]), array.array('d', [0.])
        graph.GetPoint(i, x, y)

        # Calculate the x error
        if i < n - 1:
            # Get the x value of the next point
            x_next, y_next = array.array('d', [0.]), array.array('d', [0.])
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
        graph.SetPointError(i, x_error, x_error, graph.GetErrorYlow(i), graph.GetErrorYhigh(i))

    return graph

def make_histogram_from_graph(graph):
    # Create a histogram
    hist = ROOT.TH1F("hist", "Histogram from Graph", 100, 0, 1)

    # Loop over the points in the graph
    for i in range(graph.GetN()):
        # Get the x and y values of the current point
        x, y = array.array('d', [0.]), array.array('d', [0.])
        graph.GetPoint(i, x, y)

        # Fill the histogram with the y value
        hist.Fill(y)

    return hist

def scale_correction(graph, scale, fudge=1.):
    fin = None
    hc = None
    if scale == '1015':
        fin = ROOT.TFile("eec_ratio_ch.root")
        hc = fin.Get("h1_1015")

    if scale == '1530':
        fin = ROOT.TFile("eec_ratio_ch.root")
        hc = fin.Get("h1_1530")

    # Get the number of points in the graph
    n = graph.GetN()

    # Loop over the points in the graph
    for i in range(n):
        # Get the x and y values of the current point
        x, y = array.array('d', [0.]), array.array('d', [0.])
        graph.GetPoint(i, x, y)

        scale = hc.GetBinContent(hc.FindBin(x[0]))
        # Scale the y value
        if scale <= 0:
            scale = 1.
        y_scaled = y[0] * scale * fudge

        # Set the y value
        graph.SetPoint(i, x[0], y_scaled)
        
        yel, yeh = graph.GetErrorYlow(i), graph.GetErrorYhigh(i)
        xel, xeh = graph.GetErrorXlow(i), graph.GetErrorXhigh(i)

        graph.SetPointError(i, xel, xeh, yel * fudge, yeh * fudge)

    return graph

import os  
# dname = "/Users/ploskon/Library/CloudStorage/Dropbox/2024/EECs/fromKyle/"
dname = os.path.dirname(os.path.abspath(__file__))
fname_10_15 = os.path.join(dname, "fromKyle/10to15.txt")
fname_15_30 = os.path.join(dname, "fromKyle/15to30.txt")
graph10_15 = add_x_errors(create_graph(fname_10_15))
graph15_30 = add_x_errors(create_graph(fname_15_30))
canvas = ROOT.TCanvas("canvas", "My Canvas", 800, 600)
graph10_15.Draw("AP")
graph15_30.SetMarkerColor(ROOT.kRed)
graph15_30.SetLineColor(ROOT.kRed)
graph15_30.Draw("P SAME")
canvas.Draw()

fout = ROOT.TFile("EEC_DKyle.root", "RECREATE")
graph10_15.Write("graph10_15")
gr = scale_correction(graph10_15, '1015', 0.6)
fout.cd()
gr.Write("graph10_15_scaled")
graph15_30.Write("graph15_30")
gr = scale_correction(graph15_30, '1530', 0.6)
fout.cd()
gr.Write("graph15_30_scaled")
fout.Close()


