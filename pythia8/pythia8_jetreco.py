#!/usr/bin/env python3

from __future__ import print_function
import tqdm
import argparse
import os
import numpy as np
import sys
import yasp
import cppyy
import yaml
import math
import array
import ROOT
import logging
import itertools

import heppyy.util.fastjet_cppyy
import heppyy.util.pythia8_cppyy
import heppyy.util.heppyy_cppyy

from cppyy.gbl import fastjet as fj
from cppyy.gbl import Pythia8
from cppyy.gbl.std import vector

from yasp import GenericObject
from heppyy.pythia_util import configuration as pyconf
from heppyy.util.logger import Logger
log = Logger()


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

	def write_to_yaml(self, filename):
		with open(filename, 'w') as f:	
			_out_dict = {}
			for k, v in self.__dict__.items():
				if k == 'args':
					continue
				_out_dict[k] = v
			yaml.dump(_out_dict, f, default_flow_style=False)


# singleton class for storing unique analysis names
class NamesStore(object):

	__instance = None

	def __new__(cls):
		if NamesStore.__instance is None:
			NamesStore.__instance = object.__new__(cls)
			NamesStore.__instance.names = []
		return NamesStore.__instance

	def add_unique(self, name):
		if name in self.names:
			raise ValueError(f'name {name} already exists')
		self.names.append(name)
		return name

	def add(self, name):
		requested_name = name
		while name in self.names:
			name += '_'
		self.names.append(name)
		if requested_name != name:
			log.warning(f'name {requested_name} already exists, using {name} instead')
		return name


class AnalysisBase(GenericObject):

	def __init__(self, config, **kwargs):
		if isinstance(config, str):
			self.configure_from_yaml(config)
		if isinstance(config, dict):
			self.configure_from_dict(config)
		if isinstance(config, ConfigData):
			self.configure_from_dict(config.__dict__)
		super(AnalysisBase, self).__init__(**kwargs)
		if self.name is None:
			self.name = 'AnalysisBase'
		self.name = NamesStore().add_unique(self.name)
		self._log = Logger()
		self._nev = 0
		self.init()

	def open_root_file(self):
		# open a root file if needed
		if self.root_file is None:
			self.objects = []
			if self.output is None:
				self.output = 'analysis_output.root'
				self.output.replace('.root', f'_{self.name}.root')
			if self.output.endswith('.root'):
				self.output = self.output.replace('.root', f'_{self.name}.root')
			self.output = NamesStore().add_unique(self.output)
			if self.root_file is None:
				self.root_file = ROOT.TFile(self.output, 'RECREATE')
				self._log.info(f'opened root file {self.root_file.GetName()}')
		self.root_file.cd()
   
	def init(self):
		self._log.debug('AnalysisBase::init()')

	def analyze_event(self, pythia=None, hepmc=None):
		rv = True
		if pythia:
			rv = rv and self.analyze_pythia_event(pythia)
		if hepmc:
			rv = rv and self.analyze_hepmc_event(hepmc)
		self._nev += 1
		return rv

	def analyze_pythia_event(self, pythia):
		pass

	def analyze_hepmc_event(self, hepmc):
		pass

	def __del__(self):
		if self.root_file:
			self.root_file.Write()
			self.root_file.Close()
			self._log.info(f'wrote {self.root_file.GetName()}')
		else:
			self._log.debug('no root file to write')
		self.root_file = None

class JetAnalysis(AnalysisBase):
  
	def init(self):
		super(JetAnalysis, self).init()
		self.open_root_file()
		self.root_file.cd() # done already in the base class but just to be sure
		self.tn_events 	= ROOT.TNtuple(f'tn_events_{self.name}', 	'tn_events', 	'nev:xsec:ev_weight:njets:nparts')
		self.tn_jet 		= ROOT.TNtuple(f'tn_jet_{self.name}', 		'tn_jet', 		'nev:ij:pt:eta:phi:xsec:ev_weight')
		self.tn_parts 	= ROOT.TNtuple(f'tn_parts_{self.name}', 	'tn_parts', 	'nev:ij:pt:eta:phi:pid:status:xsec:ev_weight')	
		# jet finder
		self.jet_def 			= fj.JetDefinition(fj.antikt_algorithm, self.jet_R)
		self.jet_selector = fj.SelectorPtMin(self.jet_pt_min)
		if self.part_abs_eta_max:
			self.part_selector = fj.SelectorAbsEtaMax(self.part_abs_eta_max)
		if self.jet_abs_eta_max > 0:
			if self.part_abs_eta_max and self.jet_abs_eta_max > self.part_abs_eta_max - self.jet_R:
				self.jet_abs_eta_max = self.part_abs_eta_max - self.jet_R
				self._log.warning(f'jet_abs_eta_max adjusted to {self.jet_abs_eta_max}')
			self.jet_selector = self.jet_selector * fj.SelectorAbsEtaMax(self.jet_abs_eta_max)
		if self.jet_eta_max is not None:
			self.jet_selector = self.jet_selector * fj.SelectorEtaMax(self.jet_eta_max)
		if self.jet_eta_min is not None:
			self.jet_selector = self.jet_selector * fj.SelectorEtaMin(self.jet_eta_min)
		if self.jet_pt_max > 0:
			self.jet_selector = self.jet_selector * fj.SelectorPtMax(self.jet_pt_max)
		self._log.critical(f'jet definition: {self.jet_def.description()}')
		self._log.critical(f'jet selector: {self.jet_selector.description()}')
		self._log.info(f'particle selection check: self.parts_select_final={self.parts_select_final} self.parts_select_charged={self.parts_select_charged}')

	def select_final(self, select=True):
		self.parts_select_final = select

	def select_charged(self, select=True):
		self.parts_select_charged = select

	def analyze_pythia_event(self, pythia):
		self.py_parts = []
		self.parts = vector[fj.PseudoJet]()
		if self.parts_select_final:
			if self.parts_select_charged:
				self.py_parts = [(i, p, p.id(), p.status()) for i, p in enumerate(pythia.event) if p.isFinal() and p.isCharged()]
			else:
				self.py_parts = [(i, p, p.id(), p.status()) for i, p in enumerate(pythia.event) if p.isFinal()]
		else:
			self.py_parts = [(i, p, p.id(), p.status()) for i, p in enumerate(pythia.event)]
			if self._nev < 1:
				self._log.warning(f'selecting all particles final and not final state!')
		for xp in self.py_parts:
			_psj = fj.PseudoJet(xp[1].px(), xp[1].py(), xp[1].pz(), xp[1].e())
			_psj.set_user_index(xp[0])
			self.parts.push_back(_psj)
		if self.part_selector:
			self.parts = self.part_selector(self.parts)
		if len(self.parts) < 1:
			self._log.debug(f'no particles to analyze ? self.parts_select_final={self.parts_select_final} self.parts_select_charged={self.parts_select_charged}')
			self._log.debug(f'no particles to analyze len(self.parts)={len(self.parts)}')
			return False
		self.jets = fj.sorted_by_pt(self.jet_selector(self.jet_def(self.parts)))
		self.pythia_info = Pythia8.getInfo(pythia)
		self.tn_events.Fill(self._nev, self.pythia_info.sigmaGen(), self.pythia_info.weight(), len(self.jets), len(self.parts))
		for ij, j in enumerate(self.jets):
			self.tn_jet.Fill(self._nev, ij, j.pt(), j.eta(), j.phi(), self.pythia_info.sigmaGen(), self.pythia_info.weight())
		if self.write_parts:
			for p, py_p in zip(self.parts, self.py_parts):
				ijet = [i for i, j in enumerate(self.jets) if p.user_index() in [x.user_index() for x in j.constituents()]]
				if len(ijet) < 1:
					ijet = -1
				else:
					ijet = ijet[0]
				self.tn_parts.Fill(self._nev, ijet, p.pt(), p.eta(), p.phi(), py_p[2], py_p[3], self.pythia_info.sigmaGen(), self.pythia_info.weight())
		if len(self.jets) < 1:
			return False
		return True
	
class EECAnalysis(JetAnalysis):
  
	def init(self):
		super(EECAnalysis, self).init()
		self.tn_eec = ROOT.TNtuple(f'tn_eec_{self.name}', 'tn_eec', 'nev:ij:dr:pt1:pt2:eec:ptjet:xsec:ev_weight:ptcut')
		self.pt_cuts = [0, 0.15, 1, 2]

	def analyze_pythia_event(self, pythia):
		super(EECAnalysis, self).analyze_pythia_event(pythia)
		# add the EEC analysis here
		rvs = [self.analyze_eec(pythia, ptcut) for ptcut in self.pt_cuts]
		return any(rvs)
	
   
	def analyze_eec(self, pythia, ptcut):
		# Generate all pairs from parts, excluding pairs of the same element
		# self.pairs = list(itertools.combinations(parts, 2))
		# Generate all pairs from parts, including pairs of the same element
		for ij, j in enumerate(self.jets):
			_parts_cut = [p for p in j.constituents() if p.perp() > ptcut]
			_pairs = list(itertools.product(_parts_cut, repeat=2))
			log.debug(f'number of pairs: {len(_pairs)} with ptcut: {ptcut}')
			if len(_pairs) < 1:
				continue
			for first, second in _pairs:
				dr = first.delta_R(second)
				eec = first.perp() * second.perp() / pow(j.perp(), 2.)
				self.tn_eec.Fill(self._nev, ij, dr, first.perp(), second.perp(), eec, j.perp(), self.pythia_info.sigmaGen(), self.pythia_info.weight(), ptcut)
		return True


def main():
	parser = argparse.ArgumentParser(description='read hepmc and analyze eecs', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	parser.add_argument('--config', help='configure from yaml', type=str)
	parser.add_argument('--write-config', help='write config to yaml and quit', type=str, default='')
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
	parser.add_argument('--jet-abs-eta-max', help="max eta of a jet to accept", default=0.5, type=float)	
	parser.add_argument('--part-abs-eta-max', help="max eta of a particle to accept", default=1.0, type=float)	
	parser.add_argument('--jet-R', help="jet R", default=0.4, type=float)	
	parser.add_argument('--parts-select-charged', help="only charged particles", default=False, action='store_true')
	parser.add_argument('--', dest='remainder', nargs=argparse.REMAINDER)
 
	args = parser.parse_args()	

	if args.remainder:
		log.critical(f'[i] remainder: {args.remainder}')
  
	if args.nev < args.ncounts:
		args.nev = -1
		log.info(f'[w] setting nev to -1 as nev < ncounts')

	if args.nev <= 0 and args.ncounts <= 0:
		args.nev = 1000
		log.info(f'[w] setting nev to 1000')
 
	config = ConfigData(args=args)
	if args.write_config:
		config.write_to_yaml(args.write_config)
		return
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

	D0_selector = fj.SelectorAbsEtaMax(0.8) * fj.SelectorPtMin(config.D0_min_pt) * fj.SelectorPtMax(config.D0_max_pt)
	if config.D0mode:
		log.critical(f'[i] D0 selector: {D0_selector.description()}')
 
	mycfg = []
	# pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	pythia = pyconf.create_and_init_pythia_from_args(config, mycfg)
	if not pythia:
		print("[e] pythia initialization failed.")
		return

	# jet_all = JetAnalysis(config, name='jet_all')

	eec_all = EECAnalysis(config, name='eec_all')
	eec_all.select_final(True)

	eec_ch = EECAnalysis(config, name='eec_ch')
	eec_ch.select_final(True)
	eec_ch.select_charged(True)

	_stop = False 
	pbar = tqdm.tqdm(range(config.nev))
	pbar.set_description('[i] events processed')
	pbar_counts = tqdm.tqdm(range(config.ncounts))
	pbar_counts.set_description('[i]  counts accepted')
	n_count = 0
	while not _stop:
		if not pythia.next():
			continue
		rv = True
		# rv = rv and jet_all.analyze_event(pythia=pythia)
		rv = rv and eec_all.analyze_event(pythia=pythia)
		rv = rv and eec_ch.analyze_event(pythia=pythia)
		pbar.update(1)
		if rv:
			pbar_counts.update(1)
		if pbar.n >= config.nev and config.nev > 0:
			_stop = True
		if pbar_counts.n >= config.ncounts and config.ncounts > 0:
			_stop = True
	pbar.close()
	pbar_counts.close()
	pythia.stat()
	del pythia
 
if __name__ == '__main__':
	main()
