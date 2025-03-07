#!/usr/bin/env python3

from __future__ import print_function
import tqdm
import argparse
import os
import numpy as np
import sys
import yasp
import cppyy

import sys

import heppyy.util.fastjet_cppyy
import heppyy.util.pythia8_cppyy
import heppyy.util.heppyy_cppyy

from cppyy.gbl import fastjet as fj
from cppyy.gbl import Pythia8
from cppyy.gbl.std import vector

from cppyy.gbl import heppyy

from heppyy.pythia_util import configuration as pyconf

import eech
import itertools
import ROOT

# print(cppyy.gbl.__dict__)

from heppyy.util.logger import Logger
log = Logger()

from heppyy.util.mputils import logbins, filename_safe
import itertools


class EEChistogram:
	def __init__(self, name='eec', nbins=18):
		self.ncorrel = 2
		self.nbins = int(nbins)
		self.lbins = logbins(1.e-2, 1., self.nbins)
		hname = name
		self.h = ROOT.TH1F(hname, hname, self.nbins, self.lbins)
		self.h_pairs = ROOT.TH1F(hname+'_pairs', hname+'_pairs', self.nbins, self.lbins)

	def Scale(self, weight, option="width"):
		self.h.Scale(weight, option)
		# do not scale by bin width - just a number of pairs in the bin
		self.h_pairs.Scale(weight)

	def fill(self, parts, weight=1.):
			# Generate all pairs from parts, excluding pairs of the same element
			# self.pairs = list(itertools.combinations(parts, 2))
			# Generate all pairs from parts, including pairs of the same element
			self.pairs = list(itertools.product(parts, repeat=2))
			# Calculate the list of pairs
			self.eec2 = [(first.perp() * second.perp() * weight, first.delta_R(second)) for first, second in self.pairs]
			log.debug(f'number of pairs: {len(self.eec2)}')
			# Fill the histogram
			_ = [self.h.Fill(dr, eec) for eec, dr in self.eec2]
			_ = [self.h_pairs.Fill(dr) for eec, dr in self.eec2]
   
	def Write(self):
		self.h.Write()
		self.h_pairs.Write()


class EEC2file:
	def __init__(self, output_fname='eec2file_out.root', name='eec2', args=None):
		self.output_fname = output_fname.replace('.root', '_' + filename_safe(name) + '.root')
		self.fout = ROOT.TFile(self.output_fname, 'recreate')
		log.debug(f'[i] will write to {self.fout.GetName()}')
		self.fout.cd()
		self.h = {}
		self.h['counts'] = ROOT.TH1F('counts', 'counts', 100, 0., 1000.)
		self.h['counts_ew'] = ROOT.TH1F('counts_w', 'counts_w', 100, 0., 1000.)
	
		_nbins = 18
		_lbins = logbins(1.e-2, 1., _nbins)
		self.h['z_parts_norm'] = ROOT.TH1F('z_parts_norm', 'z_parts_norm', _nbins, _lbins)
		_lbins = logbins(1.e-2, 100., _nbins)
		self.h['pt_parts_norm'] = ROOT.TH1F('pt_parts_norm', 'pt_parts_norm', _nbins, _lbins)

		self.h['eec2_pt_0.0'] = EEChistogram(name='eec2_pt0')
		self.h['eec2_pt_0.15'] = EEChistogram(name='eec2_pt0.15')
		self.h['eec2_pt_1.0'] = EEChistogram(name='eec2_pt1')
		self.h['eec2_pt_2.0'] = EEChistogram(name='eec2_pt2')

		self.h['eec2_pt_0.0_ew'] = EEChistogram(name='eec2_pt0_ew')
		self.h['eec2_pt_0.15_ew'] = EEChistogram(name='eec2_pt0.15_ew')
		self.h['eec2_pt_1.0_ew'] = EEChistogram(name='eec2_pt1_ew')
		self.h['eec2_pt_2.0_ew'] = EEChistogram(name='eec2_pt2_ew')

	def fill(self, parts, pTweight, event_weight):
		self.h['counts'].Fill(pTweight)
		self.h['counts_ew'].Fill(pTweight, event_weight)
 
		# fill the z for the object (jet if jet passed)
		_ = [self.h['z_parts_norm'].Fill(p.perp() / pTweight) for p in parts]
		_ = [self.h['pt_parts_norm'].Fill(p.perp()) for p in parts]
 
		_pTweight2 = pTweight * pTweight
		_parts = parts
		self.h['eec2_pt_0.0'].fill(_parts, weight=1./_pTweight2)
		self.h['eec2_pt_0.0_ew'].fill(_parts, weight=1./_pTweight2 * event_weight)

		_parts = [p for p in parts if p.perp() > 0.15]
		self.h['eec2_pt_0.15'].fill(_parts, weight=1./_pTweight2)
		self.h['eec2_pt_0.15_ew'].fill(_parts, weight=1./_pTweight2 * event_weight)

		_parts = [p for p in parts if p.perp() > 1.0]
		self.h['eec2_pt_1.0'].fill(_parts, weight=1./_pTweight2)
		self.h['eec2_pt_1.0_ew'].fill(_parts, weight=1./_pTweight2 * event_weight)

		_parts = [p for p in parts if p.perp() > 2.0]
		self.h['eec2_pt_2.0'].fill(_parts, weight=1./_pTweight2)
		self.h['eec2_pt_2.0_ew'].fill(_parts, weight=1./_pTweight2 * event_weight)

	def __del__(self):
		fname = self.fout.GetName()
		self.fout.cd()
		norm = self.h['counts'].Integral()
		for h in self.h.values():
			if isinstance(h, EEChistogram):
				if norm > 0:
					log.debug(f'normalizing {h.h.GetName()} by {norm}')
					h.Scale(1./norm, "width")
			if isinstance(h, ROOT.TH1F):
				if 'norm' in h.GetName():
					if norm > 0:
						if 'bw' in h.GetName():
							h.Scale(1./norm, "width")
						else:
							h.Sumw2()
							h.Scale(1./norm)
			h.Write()
		self.fout.Write()
		self.fout.Close()
		log.debug(f'[i] written {fname}.')


def match_jets_partons(jets, partons):
	# match partons to jets
	for j in jets:
		log.debug(f'jet pt: {j.pt()} eta: {j.eta()}')
		for np in range(len(partons)):
			p = partons[np]
			if j.delta_R(p) < 0.4:
				log.debug(f' - matched parton {np} with pt: {p.perp()} eta: {p.eta()}')


def main():
	parser = argparse.ArgumentParser(description='analyze hybrid with fastjet on the fly', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	group = parser.add_mutually_exclusive_group(required=True)
	group.add_argument('--pythia', help="run pythia not read hybrid", default=False, action='store_true')
	group.add_argument('-i', '--input', help='hybrid event input file', default='')
	# parser.add_argument('-n', '--nev', help='number of events', default=10, type=int)
	parser.add_argument('-w', '--wake', help='include wake particles', action='store_true', default=False)
	parser.add_argument('--ignore-wake', help='run the negative recombiner but ignore the negative index particles', action='store_true', default=False)
	parser.add_argument('-v', '--verbose', help="be verbose", default=False, action='store_true')
	parser.add_argument('--jet-min-pt', help="minimum pT jet to accept", default=20., type=float)	
	parser.add_argument('--jet-max-pt', help="max pT jet to accept", default=-1, type=float)	
	parser.add_argument('--max-eta-jet', help="max eta of a jet to accept", default=2.5, type=float)	
	parser.add_argument('--jet-R', help="jet R", default=0.4, type=float)	
	parser.add_argument('-g', '--debug', help="write debug things", default=False, action='store_true')
	parser.add_argument('--ncorrel', help='max n correlator', type=int, default=2)
	parser.add_argument('-o','--output', help='root output filename', default='eec_hybrid.root', type=str)
	parser.add_argument('--log', help="write stdout to a log file - output file will be same as -o but .log instead of .root", default=False, action='store_true')

	args = parser.parse_args()

	stdout_file = None
	if args.log:
		stdout_file_name = args.output.replace('.root', '.log')
		stdout_file = open(stdout_file_name, 'w')
		sys.stdout = stdout_file
		sys.stderr = stdout_file

	# set up logging - this uses singleton Logger
	log_level = 'DEBUG' if args.debug else 'WARNING'
	log.set_level(log_level)
	if args.verbose:
		log.enable_console()
		log.set_level('INFO')
	if args.debug:
		log.set_level('DEBUG')
  
	log.critical(args)
  
	fj.ClusterSequence.print_banner()
	print()

	# use the energy recombiner
	ner_index = 1111
	ner = heppyy.NegativeEnergyRecombiner(ner_index)
	log.critical(ner.description())

	# set up our jet definition and a jet selector
	jet_R0 = args.jet_R
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	if args.wake:
		jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0, ner)
	# jet_selector = fj.SelectorPtMin(args.jet_min_pt) * fj.SelectorAbsEtaMax(args.max_eta_jet)
	jet_selector = fj.SelectorPtMin(args.jet_min_pt) * fj.SelectorPtMax(args.jet_min_pt+20) * fj.SelectorAbsEtaMax(1 - jet_R0 * 1.05)
	if args.jet_max_pt > 0:
		jet_selector = fj.SelectorPtMin(args.jet_min_pt) * fj.SelectorPtMax(args.jet_max_pt) * fj.SelectorAbsEtaMax(1 - jet_R0 * 1.05)
	log.critical(jet_def.description())
	log.critical(jet_selector.description())

	pythia = None
	input = None
	if args.pythia:
		pythia = pyconf.create_and_init_pythia_from_args(args)
		args.output = args.output.replace('hybrid', 'pythia8')
		if 'pythia8' not in args.output:
			args.output = args.output.replace('.root', '_pythia8.root')
		if not pythia:
			log.error('pythia initialization failed.')
			return
	else:
		input = heppyy.HybridFile(args.input, medium_offset = ner_index)

	if args.nev < 10:
		args.nev = 10

	halt = EEC2file(output_fname=args.output, name='eec2', args=args)
	halt_nw	= None
	halt_ww	= None
	if args.wake:
		halt_nw = EEC2file(output_fname=args.output, name='eec2_nw', args=args)
		halt_ww = EEC2file(output_fname=args.output, name='eec2_ww', args=args)

	hnw = None
	h = eech.EEChistograms(args=args)
	log.info(h)
	if args.wake and not args.ignore_wake:
		args.output = args.output.replace('.root', '_ignore_wake.root')
		hnw = eech.EEChistograms(args=args, name='ignore_wake')

	_stop = False 
	pbar = tqdm.tqdm(range(args.nev))
	njets = 0
	while not _stop:
		if pythia:
			if not pythia.next():
				break
			parts = vector[fj.PseudoJet]([fj.PseudoJet(p.px(), p.py(), p.pz(), p.e()) for p in pythia.event if p.isFinal() and p.isCharged()])
		if input:
			if not input.nextEvent():
				break
			parts = input.getParticles(include_wake=args.wake)
			log.debug(f'number of particles: {len(parts)}')
			if args.ignore_wake:
				parts = vector[fj.PseudoJet]([p for p in input.getParticles() if p.user_index() >= 0 and p.user_index() != ner_index])
		# partons = input.getPartons()
		# sparts = input.getParticlesStr()
		# spartons = input.getPartonsStr()

		# get the event information
		ev_info = None
		ev_weight = 1.
		sigmaGen = 1.
  
		if pythia:
			ev_info = Pythia8.getInfo(pythia)
			sigmaGen = ev_info.sigmaGen()
			ev_weight = ev_info.weight()

		if input:
			ev_info = input.info()
			sigmaGen = ev_info.sigmaGen()
			ev_weight = 1.

		jets = fj.sorted_by_pt(jet_selector(jet_def(parts)))
		log.debug(f'number of jets: {len(jets)}')
		njets += len(jets)
		if jets.size() > 0:
			for j in jets:
				h.fill_jet(j, j.constituents(), j.perp(), sigmaGen=sigmaGen, weight=ev_weight)
				for p in j.constituents():
					if p.user_index() < 0:
						log.debug(f'{p.user_index()} index particle: {p.user_index()} {p.perp()} {p.eta()} {p.phi()}')
					if p.user_index() >= ner_index:
						log.debug(f'{ner_index} index particle: {p.user_index()} {p.perp()} {p.eta()} {p.phi()}')
				if hnw:
					_c = vector(fj.PseudoJet)(p for p in j.constituents() if p.user_index() >= 0)
					hnw.fill_jet(j, _c, j.perp(), sigmaGen=sigmaGen, weight=ev_weight)
				_parts = j.constituents()
				_pTweight = j.perp()
				_event_weight = ev_weight * sigmaGen
				halt.fill(_parts, _pTweight, _event_weight)
				if halt_ww:
					_parts = [p for p in j.constituents() if abs(p.user_index()) > ner_index]
					halt_ww.fill(_parts, _pTweight, _event_weight)
				if halt_nw:
					_parts = [p for p in j.constituents() if p.user_index() >= 0 and p.user_index() < ner_index]
					halt_nw.fill(_parts, _pTweight, _event_weight)

		else:
			continue
		pbar.update(jets.size())
		if pbar.n >= args.nev:
			_stop = True

	print('[i] number of jets:', njets)
	if stdout_file:
		stdout_file.close()

if __name__ == '__main__':
	main()
