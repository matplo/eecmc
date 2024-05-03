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


def get_D0s(hepmc_event):
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
		_daughters_pids = [_p.pid for _p in p.end_vertex.particles_out]
		if 211 in _daughters_pids and -321 in _daughters_pids:
				for _p in p.end_vertex.particles_out:
					Ddaughters.append(_p)
				Darray.append(p)
	return Darray, Ddaughters

def get_final_except(hepmc_event, parts, charged_only=False):
	fjparts = []
	for i,p in enumerate(hepmc_event.particles):
		if p in parts:
			continue
		if charged_only and particle.Particle.from_pdgid(p.pid).charge == 0:
				continue
		if p.status == 1 and not p.end_vertex:
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


def main():
	parser = argparse.ArgumentParser(description='read hepmc and analyze eecs', prog=os.path.basename(__file__))
	parser.add_argument('-i', '--input', help='input file', default='low', type=str, required=True)
	parser.add_argument('--hepmc', help='what format 2 or 3', default=2, type=int)
	parser.add_argument('--nev', help='number of events', default=-1, type=int)
	parser.add_argument('--ncorrel', help='max n correlator', type=int, default=2)
	parser.add_argument('-o','--output', help='root output filename', default='eec.root', type=str)
	parser.add_argument('-v', '--verbose', help="be verbose", default=False, action='store_true')
	parser.add_argument('-g', '--debug', help="write debug things", default=False, action='store_true')
	parser.add_argument('--log', help="write stdout to a log file - output file will be same as -o but .log instead of .root", default=False, action='store_true')
	parser.add_argument('--jet-min-pt', help="minimum pT jet to accept", default=20., type=float)	
	parser.add_argument('--jet-max-pt', help="max pT jet to accept", default=-1, type=float)	
	parser.add_argument('--D0-min-pt', help="minimum pT D0", default=3., type=float)	
	parser.add_argument('--D0-max-pt', help="max pT D0", default=100, type=float)	
	parser.add_argument('--max-eta-jet', help="max eta of a jet to accept", default=0, type=float)	
	parser.add_argument('--jet-R', help="jet R", default=0.4, type=float)	
	parser.add_argument('--charged-only', help="only charged particles", default=False, action='store_true')

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


	if args.nev < 0:
		args.nev = hepmc_count_events.get_n_per_file(args.input, args.hepmc)[0]

	###
	# now lets read the HEPMC file and do some jet finding
	if args.hepmc == 3:
		input_hepmc = pyhepmc.io.ReaderAscii(args.input)
	if args.hepmc == 2:
		input_hepmc = pyhepmc.io.ReaderAsciiHepMC2(args.input)

	if input_hepmc.failed():
		print ("[error] unable to read from {}".format(args.input))
		sys.exit(1)

	# jet finder
	# print the banner first
	fj.ClusterSequence.print_banner()
	print()
	jet_R0 = args.jet_R
	if args.max_eta_jet == 0:
		args.max_eta_jet = 0.9 - jet_R0 * 1.05
		log.critical(f'[i] setting max eta jet to {args.max_eta_jet}')
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	# jet_selector = fj.SelectorPtMin(args.jet_min_pt) * fj.SelectorPtMax(args.jet_max_pt) * fj.SelectorAbsEtaMax(0.9 - jet_R0 * 1.05)
	jet_selector = fj.SelectorPtMin(args.jet_min_pt) * fj.SelectorPtMax(args.jet_max_pt) * fj.SelectorAbsEtaMax(args.max_eta_jet)
	D0_selector = fj.SelectorAbsEtaMax(0.8) * fj.SelectorPtMin(args.D0_min_pt) * fj.SelectorPtMax(args.D0_max_pt)

	# from FJ contrib - not clear how to use this
	# eec = fj.contrib.EnergyCorrelator(2, 1) # default is measure pt_R
	# log.debug(eec.description())
 
	h = heec.EEC2file(args.output, name='eec2', args=args)
 
	D0indexMark = 99421
	event_hepmc = pyhepmc.GenEvent()
	pbar = tqdm.tqdm(range(args.nev))
	njets = 0
	nev_count = 0
	while not input_hepmc.failed():
		ev = input_hepmc.read_event(event_hepmc)
		if input_hepmc.failed():
			break
		nev_count += 1
		fjparts = vector[fj.PseudoJet]()
		Darray, Ddaughters = get_D0s(event_hepmc)
		if len(Darray) == 0:
			continue

		log.debug(f'---- number of D0s: {len(Darray)}')
		# get the final state particles except the D0s	
		fjparts = get_final_except(event_hepmc, Ddaughters, args.charged_only)
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
				h.fill(j.constituents(), j.perp(), event_hepmc.weight() * event_hepmc.cross_section.xsec())
  
		pbar.update(nD0jets)
		if pbar.n >= args.nev:
			break

	log.critical(f'[i] number of D0 jets accepted: {pbar.n}')
	log.critical(f'[i] number of events analyzed: {nev_count}')

if __name__ == '__main__':
	main()
