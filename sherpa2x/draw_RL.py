#!/usr/bin/env python3

from __future__ import print_function
import tqdm
import argparse
import os
import numpy as np
import sys
import yasp
import sys
import array
import ROOT

def logbins(xmin, xmax, nbins):
        lspace = np.logspace(np.log10(xmin), np.log10(xmax), nbins+1)
        arr = array.array('f', lspace)
        return arr

def main():
	parser = argparse.ArgumentParser(description='pythia8 in python', prog=os.path.basename(__file__))
	parser.add_argument('-i', '--input', help='input file', default='low', type=str, required=True)
	parser.add_argument('-o','--output', help='root output filename', default='', type=str)
	parser.add_argument('--nbins', help='number of bins in RL histogram', default=40, type=int)
	parser.add_argument('--ncorrel', help='max n correlator', type=int, default=2)
	parser.add_argument('--use-h', help='syntax f:h - use a histogram h from file f for binning', default='', type=str)
	args = parser.parse_args()

	nbins = int(args.nbins)
	lbins = logbins(1.e-3, 1., nbins)

	if args.use_h:
		try:
			_fname = args.use_h.split(':')[0]
			_hname = args.use_h.split(':')[1]
			_lbins = []
			with ROOT.TFile(_fname) as _f:
				htemplate = _f.Get(_hname)
				nbins = htemplate.GetNbinsX()
				for b in range(1, nbins+1):
					_lbins.append(htemplate.GetBinLowEdge(b))
				_lbins.append(htemplate.GetBinLowEdge(nbins) + htemplate.GetBinWidth(b))
				lbins = array.array('f', _lbins)
			print(f'[i] adopting binning {nbins} {lbins}')
		except:
			print(f'[e] unable to get the histogram {args.use_h}')
			return -1

	if not args.output:
		args.output = 'h_{}'.format(args.input)

	if args.ncorrel < 2:
		args.ncorrel = 2

	fin = ROOT.TFile(args.input)
	tn = fin.Get('tn')
	tnjet = fin.Get('tnjet')
	if not tn:
		return

	fout = ROOT.TFile(args.output, 'recreate')

	ptcuts = [0., 1.]
	# ptjetranges = [ [20, 40], [40, 60], [60, 80]]
	ptjetranges = [ [20, 40], [40, 60], [60, 80], [80, 100]]
	for nc in range(2, args.ncorrel + 1):
		for ptcut in ptcuts:
			for ptjet in ptjetranges:
				print('[i]', 'jet pt', ptjet, 'ncorrel', nc, 'ptcut', ptcut)
				fout.cd()
				hname = f'heec_{nc}_{ptjet[0]}_{ptjet[1]}_{ptcut}'.replace('.', 'p')
				cond = f'w*(jetpt > {ptjet[0]} && jetpt < {ptjet[1]} && ptpartcut=={ptcut} && n=={nc})'
				print(' drawing with condition is', cond)
				h = ROOT.TH1F(hname, hname, nbins, lbins)
				tn.Draw(f"RL>>{hname}", cond, 'e')
				hname = hname + '_jet'
				hjet = ROOT.TH1F(hname, hname, 100, 0, 100)
				cond = f'pt > {ptjet[0]} && pt < {ptjet[1]}'
				print(' drawing with condition is', cond)
				tnjet.Draw(f"pt>>{hname}", cond, 'e')
				h.Sumw2()
				njets = hjet.GetEntries()
				print(f' scale njets = {njets}')
				if njets > 0:
					h.Scale(1./njets, "width")
				else:
					h.Reset()
				fout.Write()
	fout.Close()
	print('[i] written', fout.GetName())

if __name__ == "__main__":
	main()
