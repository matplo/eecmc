#!/usr/bin/env python

import ROOT

def main():
  # make a ratio of two histograms
  f1 = ROOT.TFile("h_eec_1015_lead5.root")
  h1_1015 = f1.Get("h_eec_eecs_1015_lead5")
  h2_1015 = f1.Get("h_eec_ch_eecs_1015_lead5")
  h1_1015.Divide(h2_1015)
  
  f2 = ROOT.TFile("h_eec_1530_lead5.root")
  h1_1530 = f2.Get("h_eec_eecs_1530_lead5")
  h2_1530 = f2.Get("h_eec_ch_eecs_1530_lead5")
  h1_1530.Divide(h2_1530)

  fout = ROOT.TFile("eec_ratio_ch.root", "RECREATE")
  h1_1015.Write('h1_1015')
  h1_1530.Write('h1_1530')
  fout.Close()
  
if __name__ == '__main__':
  main()