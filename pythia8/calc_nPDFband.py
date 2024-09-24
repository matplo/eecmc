#!/usr/bin/env python3

import argparse
import os
import fnmatch
import sys
import ROOT
import numpy as np

sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, 'analysis'))
from analysis import SingleRootFile


def histogram2numpy(h):
    n = h.GetNbinsX()
    x = np.zeros(n)
    y = np.zeros(n)
    ex = np.zeros(n)
    ey = np.zeros(n)
    for i in range(n):
        x[i] = h.GetBinCenter(i)
        y[i] = h.GetBinContent(i)
        ex[i] = h.GetBinError(i)
        ey[i] = h.GetBinError(i)
    return x, y, ex, ey

def numpy2histogram(x, y, ex, ey):
    n = len(x)
    h = ROOT.TH1D('h', 'h', n, x[0], x[-1])
    for i in range(n):
        h.SetBinContent(i, y[i])
        h.SetBinError(i, ey[i])
    return h

def hist_diff_numpy(h1, h2):
    x1, y1, ex1, ey1 = histogram2numpy(h1)
    x2, y2, ex2, ey2 = histogram2numpy(h2)
    y = y1 - y2
    ey = np.sqrt(ey1**2 + ey2**2)
    return x1, y, ex1, ey

def hist_diff_signed(h1, h2):
    h_diff_plus = h1.Clone(f'{h1.GetName()}__diff_plus')
    h_diff_plus.SetDirectory(0)
    h_diff_minus = h1.Clone(f'{h1.GetName()}__diff_minus')
    h_diff_minus.SetDirectory(0)
    h_diff_plus.Add(h2, -1)
    h_diff_minus.Add(h2, -1)
    for ib in range(h_diff_plus.GetNbinsX()):
        diff = h_diff_plus.GetBinContent(ib)
        if diff < 0:
            h_diff_plus.SetBinContent(ib, 0)
        else:
            h_diff_minus.SetBinContent(ib, 0)
    return h_diff_plus, h_diff_minus
    

def find_files(rootdir='.', pattern='*'):
    return [os.path.join(rootdir, filename)
            for rootdir, dirnames, filenames in os.walk(rootdir)
            for filename in filenames
            if fnmatch.fnmatch(filename, pattern)]
    

def main():
    parser = argparse.ArgumentParser(description='process list of files with different nPDF members and produce the uncertainties', prog=os.path.basename(__file__))
    parser.add_argument('input_dir', help='input directory with root files')
    parser.add_argument('output_file', help='output directory')
    parser.add_argument('--filter', help='filter for the root files', default='*_h.root')
    parser.add_argument('--nbins', help='flag for small or std nbins', type=int, default=0) # 1 is smallnbin

    args = parser.parse_args()
    print(args)

    files = []
    _files = find_files(args.input_dir, args.filter)
    for fn in _files:
        if args.nbins == 0:
            if '_smallnbin_' in fn:
                continue
            files.append(fn)
        else:
            if '_smallnbin_' in fn:
                files.append(fn)
            
    print(files)
    print('number of files', len(files))
 
    # get the file with _0_h.root ending
    f0arr = [f for f in files if '_0_' in f]
    f0 = f0arr[0]
    print(f0arr)
    print(f0)
 
    # For each h0 iterate on the other files and get the histograms
    # create a new histogram with the sum of the squares of the differences
    # between the histograms
    # then take the square root of the sum and divide by the number of histograms
    # to get the uncertainty
    rf = SingleRootFile(args.output_file)
    rf.root_file.cd()

    # Get the list of histograms from f0
    tf0 = ROOT.TFile(f0)
    for key in tf0.GetListOfKeys():
        h0 = key.ReadObj()
        if not h0.InheritsFrom("TH1"):  # Check if the object is a histogram
            continue
        hname = h0.GetName()
        if 'h_nev_xsec_eec' in hname:
            continue
        print(hname)
        # h0.Scale(1/h0.Integral())
        rf.root_file.cd()
        h0_central = h0.Clone(f'{h0.GetName()}_central')
        h0_central.SetDirectory(0)
        h0_plus = h0.Clone(f'{h0.GetName()}_band_plus')
        h0_plus.SetDirectory(0)
        h0_minus = h0.Clone(f'{h0.GetName()}_band_minus')
        h0_minus.SetDirectory(0)
        h_diff_plus_list = []
        h_diff_minus_list = []
        for f in files:
            if f == f0:
                continue
            tf = ROOT.TFile(f)
            h = tf.Get(hname)
            if not h:
                print(f'histogram {hname} not found in {f}')
                continue
            h_diff_plus, h_diff_minus = hist_diff_signed(h, h0)
            h_diff_plus_list.append(h_diff_plus)
            h_diff_minus_list.append(h_diff_minus)
            tf.Close()
        for ib in range(h0.GetNbinsX()):
            diffs_plus = [h.GetBinContent(ib) for h in h_diff_plus_list]
            diffs_minus = [h.GetBinContent(ib) for h in h_diff_minus_list]
            h0_plus.SetBinContent(ib, np.sqrt(np.sum(np.array(diffs_plus)**2)))
            h0_minus.SetBinContent(ib, np.sqrt(np.sum(np.array(diffs_minus)**2)))
        rf.add(h0_central)
        rf.add(h0_plus)
        rf.add(h0_minus)
        # g0 = ROOT.TGraphAsymmErrors(h0)
        g0 = ROOT.TGraphAsymmErrors(h0.GetNbinsX())
        g0.SetName(f'{h0.GetName()}_band_graph')
        for ib in range(1, h0.GetNbinsX() + 1):
            g0.SetPoint(ib, h0.GetBinCenter(ib), h0.GetBinContent(ib))
            g0.SetPointEXlow(ib, h0.GetBinWidth(ib) / 2)
            g0.SetPointEXhigh(ib, h0.GetBinWidth(ib) / 2)
            g0.SetPointEYhigh(ib, h0_plus.GetBinContent(ib))
            g0.SetPointEYlow(ib, h0_minus.GetBinContent(ib))
        rf.add(g0)

        #h0.Scale(1/len(files))
        #h0.Write()

    rf.close()

if __name__ == "__main__":
    main()