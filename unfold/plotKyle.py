#!/usr/bin/env python

import ROOT
import numpy as np
import array


def logbins(xmin, xmax, nbins):
    xmin = max(xmin, 1e-10)
    xmax = max(xmax, 1e-10)
    lspace = np.logspace(np.log10(xmin), np.log10(xmax), nbins+1)
    arr = array.array('f', lspace)
    return arr


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


def find_y(x, xarr, yarr):
    min_delta = None
    idx_x_selected = -1
    for i, _x in enumerate(xarr):
        _delta = abs(_x - x)
        if min_delta is None or min_delta > _delta:
            min_delta = _delta
            idx_x_selected = i
    return yarr[idx_x_selected]
                
    
def make_histogram_from_graph_extrapolate(graph, hname, nbins = -1, logx=False):
    # Create a histogram
    n = int(graph.GetN())
    x = []
    y = []
    # Loop over the points in the graph
    for i in range(n):
        # Get the x and y values of the current point
        _x, _y = array.array('d', [0.]), array.array('d', [0.])
        graph.GetPoint(i, _x, _y)
        x.append(_x[0])
        y.append(_y[0])

    xmin = x[0]
    # xmax = x[-1] + abs(x[-1] - x[0]) * 1.01
    xmax = 1.
    print (xmin, xmax)
    if nbins < 0:
        nbins = n
    if logx:
        lbins = logbins(xmin, xmax, nbins)
        hist = ROOT.TH1F(hname, "Histogram from Graph", nbins, lbins)
    else:  
        hist = ROOT.TH1F(hname, "Histogram from Graph", nbins, xmin, xmax)

    # Loop over the points in the graph
    for i in range(nbins):
        _x = hist.GetBinCenter(i+1)
        _y = find_y(_x, x, y)
        hist.SetBinContent(i+1, _y)

    return hist

import os  
# dname = "/Users/ploskon/Library/CloudStorage/Dropbox/2024/EECs/fromKyle/"
dname = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../pQCD/fromKyle/')
fname_10_15 = os.path.join(dname, "10to15.txt")
fname_15_30 = os.path.join(dname, "15to30.txt")
graph10_15 = add_x_errors(create_graph(fname_10_15))
graph15_30 = add_x_errors(create_graph(fname_15_30))
canvas = ROOT.TCanvas("canvas", "My Canvas", 800, 600)
graph10_15.Draw("AP")
graph15_30.SetMarkerColor(ROOT.kRed)
graph15_30.SetLineColor(ROOT.kRed)
graph15_30.Draw("P SAME")
canvas.Draw()

nbins = 100
fout = ROOT.TFile("EEC_DKyle.root", "RECREATE")
graph10_15.Write("graph10_15")
h10_15 = make_histogram_from_graph_extrapolate(graph10_15, 'h_10_15', nbins)
h10_15.Write()
h10_15logx = make_histogram_from_graph_extrapolate(graph10_15, 'h_10_15_logx', nbins, True)
h10_15logx.Write()
graph15_30.Write("graph15_30")
h15_30 = make_histogram_from_graph_extrapolate(graph15_30, 'h_15_30', nbins)
h15_30.Write()
h15_30logx = make_histogram_from_graph_extrapolate(graph15_30, 'h_15_30_logx', nbins, True)
h15_30logx.Write()
fout.Close()

