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

from yasp import GenericObject


def logbins(xmin, xmax, nbins):
        lspace = np.logspace(np.log10(xmin), np.log10(xmax), nbins+1)
        arr = array.array('f', lspace)
        return arr


def find_jets_pythia(jet_def, jet_selector, pythia):
	fjparts = []
	fjparts = vector[fj.PseudoJet]([fj.PseudoJet(p.px(), p.py(), p.pz(), p.e()) for p in pythia.event if p.isFinal()])
	fjparts = vector[fj.PseudoJet](fjparts)
	print('- number of particles in the event:', fjparts.size())
	print('  leading particle pT:', fj.sorted_by_pt(fjparts)[0].perp())
	jets = jet_selector(jet_def(fjparts))
	return jets


def get_D0s(pythia, skipDstar = True):
	D0s = []
	D0daughters = []
	for i in range(pythia.event.size()):
		if abs(pythia.event[i].id()) == 421:
			Dstar = False
			# check if D* is a mother
			for j in range(pythia.event.size()):
				if pythia.event[j].mother1() == i or pythia.event[j].mother2() == i:
					if abs(pythia.event[j].id()) == 413:
						Dstar = True
						break
			if Dstar and skipDstar:
				continue
			D0s.append(i)
			for j in range(pythia.event.size()):
				if pythia.event[j].mother1() == i or pythia.event[j].mother2() == i:
					if abs(pythia.event[j].id()) == 211 or abs(pythia.event[j].id()) == 321:
						D0daughters.append(j)
	log.debug(f'D0s: {D0s}')
	log.debug(f'D0 daughters: {D0daughters}')
	return D0s, D0daughters


class ConfigData(GenericObject):
	def __init__(self, **kwargs):
		super(ConfigData, self).__init__(**kwargs)
		# configure from the arparse (including defaults)
		if self.args:
			self.configure_from_dict(self.args.__dict__)
		# configure from the scpeficied yaml file
		if self.config is not None:
			self.configure_from_yaml(self.config)		
		# configure from the command line - override values that are actually specified
		if self.args:
			for s in sys.argv:
				if s.startswith('--'):
					s = s[2:]
				for k, v in self.args.__dict__.items():
					if s != k:
						continue
					self.__setattr__(k, v)
		self.verbose = self.debug


def main():
	parser = argparse.ArgumentParser(description='read hepmc and analyze eecs', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	parser.add_argument('--config', help='configure from yaml', type=str)
	parser.add_argument('-o','--output', help='root output filename', default='eec_pythia8.root', type=str)
	parser.add_argument('--D0mode', help='set D0 mode', default=0, type=int)
	parser.add_argument('--ncounts', help='number of D0-jet counts', default=-1, type=int)
	parser.add_argument('--verbose', help="be verbose", default=False, action='store_true')
	parser.add_argument('--debug', help="write debug things", default=False, action='store_true')
	parser.add_argument('--log', help="write stdout to a log file - output file will be same as -o but .log instead of .root", default=False, action='store_true')
	parser.add_argument('--jet-min-pt', help="minimum pT jet to accept", default=20., type=float)	
	parser.add_argument('--jet-max-pt', help="max pT jet to accept", default=-1, type=float)	
	parser.add_argument('--D0-min-pt', help="minimum pT D0", default=3., type=float)	
	parser.add_argument('--D0-max-pt', help="max pT D0", default=100, type=float)	
	parser.add_argument('--max-eta-jet', help="max eta of a jet to accept", default=0.5, type=float)	
	parser.add_argument('--jet-R', help="jet R", default=0.4, type=float)	
	parser.add_argument('--charged-only', help="only charged particles", default=False, action='store_true')

	args = parser.parse_args()	
 
	config = ConfigData(args=args)
	log.critical(f'config: {config}')

	if args.debug:
		log.set_level(logging.DEBUG)

	if args.output == 'eec_pythia8.root':
		if args.py_vincia:
			args.output = args.output.replace('.root', '_vincia.root')
		if args.py_dire:
			args.output = args.output.replace('.root', '_dire.root')
		print("[w] using [modified] default output file:", args.output)
	else:
		print("[w] using specified output file:", args.output)

	pythia = Pythia8.Pythia()

	# jet finder
	# print the banner first
	fj.ClusterSequence.print_banner()
	print()
	jet_R0 = 0.4
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	jet_selector = fj.SelectorPtMin(args.py_pthatmin)
	jet_selector = fj.SelectorPtMin(args.py_pthatmin) * fj.SelectorPtMax(args.py_pthatmin+20) * fj.SelectorAbsEtaMax(1 - jet_R0 * 1.05)
 
	mycfg = []
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	if not pythia:
		print("[e] pythia initialization failed.")
		return

	log.critical(f'[i] jet definition: {jet_def.description()}')
	log.critical(f'[i] jet selector: {jet_selector.description()}')
	log.critical(f'[i] D0 selector: {D0_selector.description()}')

	_stop = False 
	pbar = tqdm.tqdm(range(args.nev))
	njets = 0
	while not _stop:
		if not pythia.next():
			continue
		if args.D0mode:
			D0s, D0daughters = get_D0s(pythia)
			if len(D0s) == 0:
				continue
		fjparts = vector[fj.PseudoJet]([fj.PseudoJet(p.px(), p.py(), p.pz(), p.e()) for p in pythia.event if p.isFinal() and p.isCharged()])
		# jets = fj.sorted_by_pt(jet_def(fjparts))
		jets = fj.sorted_by_pt(jet_selector(jet_def(fjparts)))
		njets += len(jets)
		# _info = Pythia8.pythia.info
		_info = Pythia8.getInfo(pythia)
		sigmaGen = _info.sigmaGen()
		ev_weight = _info.weight()
		if jets.size() > 0:
			for j in jets:
				pass
		else:
			continue
		pbar.update(jets.size())
		if pbar.n >= args.nev:
			_stop = True

	print('[i] number of jets:', njets)

if __name__ == '__main__':
	main()
