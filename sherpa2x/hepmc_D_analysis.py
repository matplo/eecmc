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
import heppyy.util.heppyy_cppyy

from cppyy.gbl import fastjet as fj
from cppyy.gbl.std import vector
from cppyy.gbl import EnergyCorrelators

import pyhepmc
import particle
import ROOT
import math
import array
import heec
import hepmc_count_events
from yasp import GenericObject

from heppyy.util.logger import Logger
log = Logger()

def logbins(xmin, xmax, nbins):
		lspace = np.logspace(np.log10(xmin), np.log10(xmax), nbins+1)
		arr = array.array('f', lspace)
		return arr


def find_jets_hepmc(jet_def, jet_selector, hepmc_event):
	fjparts = []
	# parts = vector[fj.PseudoJet]([fj.PseudoJet(p.px(), p.py(), p.pz(), p.e()) for p in pythia.event if p.isFinal()])
	for i,p in enumerate(hepmc_event.particles):
		if p.status == 1 and not p.end_vertex:
			psj = fj.PseudoJet(p.momentum.px, p.momentum.py, p.momentum.pz, p.momentum.e)
			# psj.set_user_index(i)
			fjparts.append(psj)
	fjparts = vector[fj.PseudoJet](fjparts)
	log.debug('- number of particles in the event:', fjparts.size())
	log.debug('  leading particle pT:', fj.sorted_by_pt(fjparts)[0].perp())
	jets = jet_selector(jet_def(fjparts))
	return jets


def get_D0s_skipDstar(hepmc_event, noDstar=False):
		Darray = []
		Ddaughters = []
		for i,p in enumerate(hepmc_event.particles):
				if abs(p.pid) != 421:
						continue
				if len(p.end_vertex.particles_out) != 2:
						continue
				_daughters_status = [_p.status for _p in p.end_vertex.particles_out]
				if [1, 1] != _daughters_status:
						continue
				log.debug(f'  - D0 found: {p}')
				log.debug(f'    - mothers: {p.end_vertex.particles_in}')
				log.debug(f'    - daughters: {p.end_vertex.particles_out}')
				_daughters_pids = [_p.pid for _p in p.end_vertex.particles_out]
				# Check if D0 is a result of a D* decay
				_mothers_pids = [_p.pid for _p in p.production_vertex.particles_in]
				if 413 in _mothers_pids or -413 in _mothers_pids:
						log.debug(f'  - D0 is a result of a D* decay - skipping')
						for _p in p.production_vertex.particles_in:
							log.debug(f' 	- D* daughter: {_p} {_p.pid} {_p.status}')
						continue
				if (_daughters_pids[0] == 211 and _daughters_pids[1] == -321) or (_daughters_pids[0] == -321 and _daughters_pids[1] == 211):
						for _p in p.end_vertex.particles_out:
								Ddaughters.append(_p)
						Darray.append(p)
		return Darray, Ddaughters
	

def get_D0s(hepmc_event, noDstar=False):
		Darray = []
		Ddaughters = []
		for i, p in enumerate(hepmc_event.particles):
			# Check if D*
			if abs(p.pid) == 413:
				log.debug(f'  - D* found {p.pid} {p.status}')
				_Dstar_daughters_pids = [_p.pid for _p in p.end_vertex.particles_out]
				# if D0 decay is D* -> D0 pi
				if 421 == abs(_Dstar_daughters_pids[0]) or 421 == abs(_Dstar_daughters_pids[1]):
					# skip the pion
					if 421 == abs(_Dstar_daughters_pids[0]):
						Ddaughters.append(p.end_vertex.particles_out[1])
					if 421 == abs(_Dstar_daughters_pids[1]):
						Ddaughters.append(p.end_vertex.particles_out[0])
					for _p in p.end_vertex.particles_out:
						log.debug(f' 	- D* daughter: {_p} {_p.pid} {_p.status}')
					log.debug(f' 	- D* daughter to be skipped: {Ddaughters[-1].pid} {Ddaughters[-1].status}')
				continue
			if abs(p.pid) != 421:
				continue
			if len(p.end_vertex.particles_out) != 2:
				continue
			_daughters_status = [_p.status for _p in p.end_vertex.particles_out]
			if [1, 1] != _daughters_status:
				continue
			log.debug(f'  - D0 found: {p}')
			log.debug(f'    - mothers: {p.end_vertex.particles_in}')
			log.debug(f'    - daughters: {p.end_vertex.particles_out}')
			_daughters_pids = [_p.pid for _p in p.end_vertex.particles_out]
			if (_daughters_pids[0] == 211 and _daughters_pids[1] == -321) or (_daughters_pids[0] == -321 and _daughters_pids[1] == 211):
				for _p in p.end_vertex.particles_out:
					Ddaughters.append(_p)
				Darray.append(p)
		return Darray, Ddaughters


def get_final_except(hepmc_event, parts, charged_only=False):
	fjparts = []
	for i,p in enumerate(hepmc_event.particles):
		if p in parts:
			continue
		if p.status == 1 and not p.end_vertex:
			if charged_only:
				_charge = None
				try: 
					_charge = particle.Particle.from_pdgid(p.pid).charge
				except particle.particle.particle.ParticleNotFound as e:
					print('[w] particle not found:', e)
					continue
				if _charge == 0:
					continue
			psj = fj.PseudoJet(p.momentum.px, p.momentum.py, p.momentum.pz, p.momentum.e)
			psj.set_user_index(i)
			fjparts.append(psj)
	fjparts = vector[fj.PseudoJet](fjparts)
	return fjparts

def print_debug(Darray, Ddaughters):
	for _p in Darray:
		log.debug(f'D0 is: {_p}')
		for _pp in _p.end_vertex.particles_out:
			log.debug(f'          {_pp}')
		for _pp in Ddaughters:
			log.debug(f'  - check:{_pp}')


class ConfigData(GenericObject):
	def __init__(self, **kwargs):
		super(ConfigData, self).__init__(**kwargs)
		if self.args:
			self.configure_from_dict(self.args.__dict__)
		self.verbose = self.debug

class Config(yasp.GenericObject):
	jet_min_pt = 20.
	jet_max_pt = -1
	D0_min_pt = 3.
	D0_max_pt = 100.
	max_eta_jet = 0
	jet_R = 0.4
	charged_only = False
	ncorrel = 2
	def __init__(self, **kwargs):
		super(Config, self).__init__(**kwargs)

def process_file(input_file, output_file, config, args):
	pass
	

def main():
	parser = argparse.ArgumentParser(description='read hepmc and analyze eecs', prog=os.path.basename(__file__))
	parser.add_argument('--input', help='input file - hepmc', type=str)
	parser.add_argument('--config', help='configure from yaml', type=str)
	parser.add_argument('--hepmc', help='what format 2 or 3', default=2, type=int)
	parser.add_argument('--nev', help='number of events', default=-1, type=int)
	parser.add_argument('--ncounts', help='number of D0-jet counts', default=-1, type=int)
	parser.add_argument('--ncorrel', help='max n correlator', type=int, default=2)
	parser.add_argument('--output', help='root output filename', default='eec.root', type=str)
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
	parser.add_argument('--use-h', help='syntax f:h - use a histogram h from file f for binning', default='', type=str)

	args = parser.parse_args()	
    # Check that at least one of --input or --config is provided
	if args.input is None and args.config is None:
		parser.error('One of --input or --config arguments is required')
        
	log.critical(args)

	config = ConfigData(args=args)
	if args.config is not None:
		config.configure_from_yaml(args.config)		
		log.critical(config)
		for s in sys.argv:
			if s.startswith('--'):
				s = s[2:]
			for k, v in args.__dict__.items():
				if s != k:
					continue
				config.__setattr__(k, v)
	log.critical(config)
 
	if config.input is None:
		log.critical('no input file provided')
		return 1

	stdout_file = None
	if config.log:
		stdout_file_name = config.output.replace('.root', '.log')
		stdout_file = open(stdout_file_name, 'w')
		sys.stdout = stdout_file
		sys.stderr = stdout_file

	# set up logging - this uses singleton Logger
	log_level = 'DEBUG' if config.debug else 'WARNING'
	log.set_level(log_level)
	if config.verbose:
		log.enable_console()
		log.set_level('INFO')
	if config.debug:
		log.set_level('DEBUG')  

	if config.nev < 0:
		config.nev = hepmc_count_events.get_n_per_file(config.input, config.hepmc)[0]

	###
	# now lets read the HEPMC file and do some jet finding
	if config.hepmc == 3:
		input_hepmc = pyhepmc.io.ReaderAscii(config.input)
	if config.hepmc == 2:
		input_hepmc = pyhepmc.io.ReaderAsciiHepMC2(config.input)

	if input_hepmc.failed():
		print ("[error] unable to read from {}".format(config.input))
		sys.exit(1)

	# jet finder
	# print the banner first
	fj.ClusterSequence.print_banner()
	print()
	jet_R0 = config.jet_R
	if config.max_eta_jet == 0:
		config.max_eta_jet = 0.9 - jet_R0 * 1.05
		log.critical(f'[i] setting max eta jet to {config.max_eta_jet}')
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	# jet_selector = fj.SelectorPtMin(config.jet_min_pt) * fj.SelectorPtMax(config.jet_max_pt) * fj.SelectorAbsEtaMax(0.9 - jet_R0 * 1.05)
	jet_selector = fj.SelectorPtMin(config.jet_min_pt) * fj.SelectorPtMax(config.jet_max_pt) * fj.SelectorAbsEtaMax(config.max_eta_jet)
	D0_selector = fj.SelectorAbsEtaMax(0.8) * fj.SelectorPtMin(config.D0_min_pt) * fj.SelectorPtMax(config.D0_max_pt)

	log.critical(f'[i] jet definition: {jet_def.description()}')
	log.critical(f'[i] jet selector: {jet_selector.description()}')
	log.critical(f'[i] D0 selector: {D0_selector.description()}')

	h = heec.EEC2file(config.output, name='eec2', load=False, use_h_binning=config.use_h, args=args)
 
	D0indexMark = 99421
	event_hepmc = pyhepmc.GenEvent()
	pbar = tqdm.tqdm(range(config.nev))
	pbar.set_description('[i] events processed')
	pbarD0 = tqdm.tqdm(range(config.ncounts))
	pbarD0.set_description('[i]   D-jet accepted')
	njets = 0
	nev_count = 0
	nD0jets_total = 0
	while not input_hepmc.failed():
		ev = input_hepmc.read_event(event_hepmc)
		if input_hepmc.failed():
			break
		nev_count += 1
		pbar.update(1)
		fjparts = vector[fj.PseudoJet]()
		Darray, Ddaughters = get_D0s(event_hepmc)
		if len(Darray) == 0:
			continue

		log.debug(f'---- number of D0s: {len(Darray)}')
		# get the final state particles except the D0s	
		fjparts = get_final_except(event_hepmc, Ddaughters, config.charged_only)
		if len(fjparts) == 0:
			continue
		# check if D0's accepted
		D0accept = False
		for p in Darray:
			psjD0 = fj.PseudoJet(p.momentum.px, p.momentum.py, p.momentum.pz, p.momentum.e)
			if not D0_selector(psjD0):
				continue
			D0accept = True
			psjD0.set_user_index(D0indexMark)
			log.debug(f'  - D0 accepted: pT={psjD0.perp()} eta={psjD0.eta()}')
			fjparts.push_back(psjD0)
	 
		if not D0accept:
			continue

		jets = fj.sorted_by_pt(jet_selector(jet_def(fjparts)))
		if len(jets) == 0:
			continue
		njets += len(jets)
		nD0jets = 0
		for j in jets:
			for _p in j.constituents():
				if _p.user_index() != D0indexMark:
					continue
				print_debug(Darray, Ddaughters)
				log.debug(f'  - jet: {j.perp()} {j.eta()} {j.phi()}')
				log.debug(f'  	- jet constituent: {_p.perp()} {_p.eta()} {_p.phi()} {_p.user_index()}')
				nD0jets += 1
				log.debug(f' event weight: {event_hepmc.weight()} cross section: {event_hepmc.cross_section.xsec()}')
				h.fill(j.constituents(), j.perp(), event_hepmc.cross_section.xsec(), event_hepmc.weight())
	
		# pbar.update(nD0jets)
		nD0jets_total += nD0jets
		pbarD0.update(nD0jets)
		if config.ncounts > 0 and nD0jets_total >= config.ncounts:
			break

	pbar.close()
	pbarD0.close()

	log.critical(f'[i] number of D0 jets accepted: {nD0jets_total}')
	log.critical(f'[i] number of events analyzed: {nev_count}')
 
	h.Write()
	h.Close()

if __name__ == '__main__':
	main()
