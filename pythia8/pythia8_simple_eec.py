#!/usr/bin/env python3

from __future__ import print_function
import tqdm
import argparse
import os
import numpy as np
import sys
import yasp
import cppyy


import heppyy.util.fastjet_cppyy
import heppyy.util.pythia8_cppyy
import heppyy.util.heppyy_cppyy

from cppyy.gbl import fastjet as fj
from cppyy.gbl import Pythia8
from cppyy.gbl.std import vector

# from cppyy.gbl import pythiaext

from heppyy.pythia_util import configuration as pyconf
import logging
from heppyy.util.logger import Logger
log = Logger()

import ROOT
import math
import array
import itertools


sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, 'analysis'))
from analysis import ConfigData, SingleRootFile


def logbins(xmin, xmax, nbins):
        lspace = np.logspace(np.log10(xmin), np.log10(xmax), nbins+1)
        arr = array.array('f', lspace)
        return arr


def main():
	parser = argparse.ArgumentParser(description='read hepmc and analyze eecs', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	parser.add_argument('--ncorrel', help='max n correlator', type=int, default=2)
	parser.add_argument('-o','--output', help='root output filename', default='pythia8_simple_eec_output.root', type=str)
	parser.add_argument('--charged', help='charged particles only', action='store_true')
	parser.add_argument('--jet-pt-min', help='jet pT min', type=float, default=20.0)
	parser.add_argument('--jet-pt-max', help='jet pT min', type=float, default=1000.0)
	parser.add_argument('--debug', help='debug', action='store_true')
	parser.add_argument('--write-config', help='write config to yaml and quit', type=str, default='')
	args = parser.parse_args()	
 
	if args.debug:
		log.set_level(logging.DEBUG)

	config = ConfigData(args=args)
	if args.write_config:
			config.write_to_yaml(args.write_config)
			return
	if config.nev <= 0 and config.ncounts <= 0:
			config.nev = 10
			log.info(f"[w] setting nev to {config.nev}")
	log.info(f"config: {config}")

	if args.output == 'pythia8_simple_eec_output.root':
		if args.py_vincia:
			args.output = args.output.replace('.root', '_vincia.root')
		if args.py_dire:
			args.output = args.output.replace('.root', '_dire.root')
		print("[w] using [modified] default output file:", args.output)
	else:
		print("[w] using specified output file:", args.output)

	rf = SingleRootFile(fname=config.output)
	rf.root_file.cd()
	tn_events 	= ROOT.TNtuple('tn_events', 'tn_events', 'nev:xsec:ev_weight:nparts:x1:x2:QFac')
	tn_jet = ROOT.TNtuple(f'tn_jet', 'tn_jet', 'nev:xsec:ev_weight:nj:ij:pt:eta:phi:m:ptlead')
	pt_cuts = [0.15, 1.0]
	tn_eec = {}
	for ptcut in pt_cuts:
		tn_eec[ptcut] = ROOT.TNtuple(f'tn_eec_ptcut{ptcut}', 'tn_eec_ptcut{ptcut}', 'nev:xsec:ev_weight:ij:dr:pt1:pt2:eec:ptjet:ptcut:iidx:jidx:ptlead')

	pythia = Pythia8.Pythia()

	# jet finder
	# print the banner first
	fj.ClusterSequence.print_banner()
	print()
	jet_R0 = 0.4
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	jet_selector = fj.SelectorPtMin(args.jet_pt_min)
	jet_selector = fj.SelectorPtMin(args.jet_pt_min) * fj.SelectorPtMax(args.jet_pt_max) * fj.SelectorAbsEtaMax(1 - jet_R0 * 1.05)
 
	mycfg = []
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	if not pythia:
		print("[e] pythia initialization failed.")
		return

	_stop = False 
	pbar = tqdm.tqdm(range(args.nev))
	njets = 0
	iev = 0
	while not _stop:
		if not pythia.next():
			continue
		iev += 1

		if args.charged:
			fjparts = vector[fj.PseudoJet]([fj.PseudoJet(p.px(), p.py(), p.pz(), p.e()) for p in pythia.event if p.isFinal() and p.isCharged()])
		else:
			fjparts = vector[fj.PseudoJet]([fj.PseudoJet(p.px(), p.py(), p.pz(), p.e()) for p in pythia.event if p.isFinal()])

		jets = fj.sorted_by_pt(jet_selector(jet_def(fjparts)))
		njets += len(jets)
		# _info = Pythia8.pythia.info - defunct
		_info = Pythia8.getInfo(pythia)
		sigmaGen = _info.sigmaGen()
		ev_weight = _info.weight()
		if jets.size() > 0:
			# tn_events 	= ROOT.TNtuple('tn_events', 'tn_events', 'nev:xsec:ev_weight:nparts')
			tn_events.Fill(iev, sigmaGen, ev_weight, pythia.event.size(), _info.x1(), _info.x2(), _info.QFac())
			for ij, j in enumerate(jets):
				ptlead = fj.sorted_by_pt(j.constituents())[0].perp()
				# print(j.perp())
				# tn_jet = ROOT.TNtuple(f'tn_jet', 'tn_jet', 'nev:xsec:ev_weight:nj:ij:pt:eta:phi:m:ptlead')				
				tn_jet.Fill(iev, sigmaGen, ev_weight, njets, ij, j.perp(), j.eta(), j.phi(), j.m(), fj.sorted_by_pt(j.constituents())[0].perp())
				for ptcut in pt_cuts:
						_parts_cut = [p for p in j.constituents() if p.perp() >= ptcut]
						_pairs = list(itertools.product(_parts_cut, repeat=2))
						if len(_pairs) < 1:
								continue
						for first, second in _pairs:
								dr = first.delta_R(second)
								eec = first.perp() * second.perp() / pow(j.perp(), 2.)
								# tn_eec[ptcut] = ROOT.TNtuple(f'tn_eec_ptcut{ptcut}', 'tn_eec_ptcut{ptcut}', 'nev:xsec:ev_weight:ij:dr:pt1:pt2:eec:ptjet:ptcut:iidx:jidx:ptlead')
								tn_eec[ptcut].Fill(iev, sigmaGen, ev_weight, ij, dr, first.perp(), second.perp(), eec, j.perp(), ptcut, first.user_index(), second.user_index(), ptlead)
		else:
			continue
		pbar.update(jets.size())
		if pbar.n >= args.nev:
			_stop = True

	pythia.stat()
 
	print('[i] number of jets:', njets)
	rf.close()

if __name__ == '__main__':
	main()
