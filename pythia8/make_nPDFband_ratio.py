#!/usr/bin/env python3

pPb_files=[ 'pythia_std_pPb_5TeV_nPDF_EPPS21nlo_CT18Anlo_Pb208.root', 
            'pythia_argantyr_pPb_5TeV_nPDF_EPPS21nlo_CT18Anlo_Pb208.root',
            '/Users/ploskon/devel/eecmc/pythia8/simple_eec_pythia_pPb_5TeV_nPDF_output_smallnbin_h.root',
            '/Users/ploskon/devel/eecmc/pythia8/simple_eec_pythia_pPb_5TeV_nPDF_argantyr_output_smallnbin_h.root']

pp_file='simple_eec_pythia_pp_5TeV_output_smallnbin_h.root'

pPb_files=[ 'pythia_std_pPb_5TeV_nPDF_EPPS21nlo_CT18Anlo_Pb208.root', 
            'pythia_argantyr_pPb_5TeV_nPDF_EPPS21nlo_CT18Anlo_Pb208.root',
            '/Users/ploskon/devel/eecmc/pythia8/simple_eec_pythia_pPb_5TeV_nPDF_output_h.root',
            '/Users/ploskon/devel/eecmc/pythia8/simple_eec_pythia_pPb_5TeV_nPDF_argantyr_output_h.root']

pp_file='simple_eec_pythia_pp_5TeV_output_h.root'

import os
import sys
import ROOT
import ctypes

sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, 'analysis'))
from analysis import SingleRootFile

def tgrDivHist(tgr, hist, newname):
    tgr = tgr.Clone(newname)
    for i in range(tgr.GetN()):
        x = ctypes.c_double(0)
        y = ctypes.c_double(0)
        tgr.GetPoint(i, x, y)
        # bin_content = hist.GetBinContent(i + 1)
        bin_content = hist.GetBinContent(hist.FindBin(x.value))
        if bin_content != 0:
            tgr.SetPoint(i, x.value, y.value / bin_content)
            tgr.SetPointEYhigh(i, tgr.GetErrorYhigh(i) / bin_content)
            tgr.SetPointEYlow(i, tgr.GetErrorYlow(i) / bin_content)
        else:
            tgr.SetPoint(i, x.value, 0)
            tgr.SetPointEYhigh(i, 0)
            tgr.SetPointEYlow(i, 0)
            print(f"Warning: Bin content is zero for bin {i + 1} in histogram {hist.GetName()}")
    return tgr

print(pp_file)
print(pPb_files)

rout = SingleRootFile('pPb_nPDFband_ratio.root')

f_pp = ROOT.TFile(pp_file)
for key in f_pp.GetListOfKeys():
    hname = key.GetName()
    h = key.ReadObj()
    if not h.InheritsFrom('TH1'):
        continue
    if 'h_eec' not in hname:
        continue
    print('Processing', hname)
    gr_name0 = hname
    gr_name1 = hname + '_band_graph'
    for fn in pPb_files:
        # print(' -> ', fn, 'for', hname)
        f_pPb = ROOT.TFile(fn)
        if f_pPb.IsZombie():
            print(f"Error opening file {fn}")
            continue
        gr_name = gr_name0
        _k = f_pPb.Get(gr_name)
        if not _k:
            gr_name = gr_name1
            _k = f_pPb.Get(gr_name1)
            if not _k:
                f_pPb.Close()
                continue

        print(f"Processing {gr_name} in {fn}")
        if 'argantyr' in fn:
            _extraname = '_argantyr'
        else:
            _extraname = ''
        if _k.InheritsFrom('TGraphAsymmErrors'):
            rout.root_file.cd()
            gr_ratio = tgrDivHist(_k, h, hname + '_eec_nPDFband_ratio' + _extraname)
            gr_ratio.SetName(hname + '_eec_nPDFband_ratio' + _extraname)
            gr_ratio.SetTitle(fn)
            rout.add(gr_ratio)
        elif _k.InheritsFrom('TH1'):
            rout.root_file.cd()
            gr_ratio = ROOT.TGraphAsymmErrors(_k, h, 'pois')
            gr_ratio.SetName(hname + '_eec_ratio' + _extraname)
            gr_ratio.SetTitle(fn)
            rout.add(gr_ratio)
        f_pPb.Close()
        rout.root_file.cd()
       # rout.root_file.cd()
       # gr_ratio.Write()
        f_pPb.Close()
f_pp.Close()

rout.close()