#!/usr/bin/env python3

import ROOT
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, 'analysis'))
from analysis import SingleRootFile

pp_files = [
"/Users/ploskon/devel/eecmc/pythia8/pcone_bg_pp_h.root 		:h_eec0_eec :p  : title=p_{T, part} > 0.15 GeV/c -- pp",
"/Users/ploskon/devel/eecmc/pythia8/pcone_bg_pp_h.root 		:h_eec0_pconebg_eec :hist -k +f1001 +a20 : title=perp. cone background",
"/Users/ploskon/devel/eecmc/pythia8/pcone_bg_pp_h.root 		:h_eec1_eec :p  : title=p_{T, part} > 1.0 GeV/c -- pp",
"/Users/ploskon/devel/eecmc/pythia8/pcone_bg_pp_h.root 		:h_eec1_pconebg_eec :hist -k +f1001 +a20: title=perp. cone background"]

pPb_files = [
"/Users/ploskon/devel/eecmc/pythia8/pcone_bg_pPb_argantyr_h.root 		:h_eec0_eec :p  : title=p_{T, part} > 0.15 GeV/c -- pPb Argantyr",
"/Users/ploskon/devel/eecmc/pythia8/pcone_bg_pPb_argantyr_h.root 		:h_eec0_pconebg_eec :hist -k +f1001 +a20 : title=perp. cone background",
"/Users/ploskon/devel/eecmc/pythia8/pcone_bg_pPb_argantyr_h.root 		:h_eec1_eec :p  : title=p_{T, part} > 1.0 GeV/c -- pPb Argantyr",
"/Users/ploskon/devel/eecmc/pythia8/pcone_bg_pPb_argantyr_h.root 		:h_eec1_pconebg_eec :hist -k +f1001 +a20: title=perp. cone background"]

file_modif = None
file_modif = ['_h.root', '_smallnbin_h.root']
              
def get_hist(s, syst=''):
    fin = s.split(":")[0].strip()
    if file_modif:
      fin = fin.replace(file_modif[0], file_modif[1])
    hin = s.split(":")[1].strip()
    print(fin, hin)
    f = ROOT.TFile(fin)
    h = f.Get(hin)
    h.SetName(f'{h.GetName()}_{syst}')
    h.SetDirectory(0)
    f.Close()
    return h
    
def subtract_bg(h_eec, h_bg):
    h_eec_bg = h_eec.Clone(f'{h_eec.GetName()}_bgsub')
    h_eec_bg.SetDirectory(0)
    h_eec_bg.Add(h_bg, -1)
    return h_eec_bg
  
def divide(h1, h2):
    h = h1.Clone(f'{h1.GetName()}_div_{h2.GetName()}')
    h.SetDirectory(0)
    h.Divide(h2)
    return h 

def main():
    ROOT.gROOT.SetBatch()
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    
    fout = SingleRootFile("subtract_pcone_bg_ratio.root")
    fout.root_file.cd()
    
    # first subtract the background
    pp_hists = []
    for i in range(0, len(pp_files), 2):
        h_eec = get_hist(pp_files[i], 'pp')
        pp_hists.append(h_eec)        
        h_bg = get_hist(pp_files[i+1], 'pp')
        #pp_hists.append(h_bg)
        h_eec_bg = subtract_bg(h_eec, h_bg)
        pp_hists.append(h_eec_bg)        

    for h in pp_hists:
      print(h.GetName())

    pPb_hists = []
    for i in range(0, len(pPb_files), 2):
        h_eec = get_hist(pPb_files[i], 'pPb')
        pPb_hists.append(h_eec)        
        h_bg = get_hist(pPb_files[i+1], 'pPb')
        #pPb_hists.append(h_bg)
        h_eec_bg = subtract_bg(h_eec, h_bg)
        pPb_hists.append(h_eec_bg)        

    for h in pPb_hists:
      print(h.GetName())

    ratios = []        
    # then divide the pp and pPb histograms
    for i in range(4):
        h_eec_pp = pp_hists[i]
        h_eec_pPb = pPb_hists[i]
        h_ratio = divide(h_eec_pPb, h_eec_pp)
        ratios.append(h_ratio)        

    for h in ratios:
      print(h.GetName())

    for h in pp_hists:
        fout.add(h)
    for h in pPb_hists:
        fout.add(h)
    for h in ratios:
        fout.add(h)
    fout.close()
    
if __name__ == "__main__":      
    main() 