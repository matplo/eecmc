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
import pyhepmc

import heppyy.util.fastjet_cppyy
import heppyy.util.pythia8_cppyy
import heppyy.util.heppyy_cppyy

from cppyy.gbl import fastjet as fj
from cppyy.gbl import Pythia8
from cppyy.gbl.std import vector

from heppyy.pythia_util import configuration as pyconf
from heppyy.util.logger import Logger
log = Logger()

sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, 'analysis'))
sys.path.append(os.path.join(os.path.dirname(__file__), 'analysis'))
from analysis import ConfigData, JetAnalysis, EECAnalysis, EECAnalysisCorrel, JetChargedFullAnalysis, DataSource, DataSourceD0, SingleRootFile
import hepmc_count_events

def logbins(xmin, xmax, nbins):
        lspace = np.logspace(np.log10(xmin), np.log10(xmax), nbins+1)
        arr = array.array('f', lspace)
        return arr


def is_D0_daughter(p, event):
  mothers = p.motherList()
  for m in mothers:
    if abs(event[m].id()) == 421:
      return event[m]
  return False

def is_D0_Kpi_decay(p, event):
  if abs(p.id()) == 421:
    daughters = p.daughterList()
    if len(daughters) != 2:
      return False
    if abs(event[daughters[0]].id()) == 321 and abs(event[daughters[1]].id()) == 211:
      return [daughters[0], daughters[1]]
    if abs(event[daughters[0]].id()) == 211 and abs(event[daughters[1]].id()) == 321:
      return [daughters[1], daughters[0]]
  return False

def D0_any_decay(p, event):
  if abs(p.id()) == 421:
    daughters = p.daughterList()
    _ds = [i for i in daughters]
    return _ds
  return False
    
def is_Dstar_daughter(p, event):
  mothers = p.motherList()
  for m in mothers:
    if abs(event[m].id()) == 413:
      return event[m]
  return False

def is_D0_daughter(p, event):
  mothers = p.motherList()
  for m in mothers:
    if abs(event[m].id()) == 421:
      return event[m]
  return False

# D0 modes:
# 0 - no D0 selection
# 1 - D0 selection
# 2 - D0 selection and remove D0 daughters
# 3 - D0 selection Kpi decay
# 4 - D0 selection Kpi decay and remove D0 daughters

# important switch: veto D* daughters --vetoDstar

# class that operates on pythia8 events
# find jets on parton level
# find jets on hadron level using only charged visible hadrons
# match the hadron and parton jets via _deltaR(R) criterion
# calculate EEC for both parton and hadron jets
# calculate EEC for charged particles in the hadron jets
# write the results to a root file as ntuple
# the ntuple has the following structure:
# ptjet, ptjetch, ptlead, ptleadch, dr, eec
# ptjet - pT of the jet
# ptjetch - pT of the charged jet
# ptlead - pT of the leading particle in the jet
# ptleadch - pT of the leading charged particle in the jet
# dr - delta R between the particles
# eec - EEC value
# the ntuple is filled for each jet in the event
# the event is selected if at least one jet is found

class Pythia8JetEECAnalysisCorrel(object):
  def __init__(self, name, config):
    self.config = config

    self.D0_selector = fj.SelectorAbsRapMax(self.config.D0_abs_eta_max) * fj.SelectorPtMin(self.config.D0_pt_min) * fj.SelectorPtMax(self.config.D0_pt_max)
    self.ch_hadron_selector = fj.SelectorPtMin(self.config.part_pt_min) * fj.SelectorAbsEtaMax(self.config.part_abs_eta_max)

    self.rf = SingleRootFile(fname=config.output)
    self.rf.root_file.cd()
    stn = f'tn_{self.name}'
    self.tn_jet = ROOT.TNtuple(f'{stn}_jet', 'tn_jet', 'ptjet:ptjetch:ptlead:ptleadch')
    self.tn_eec = ROOT.TNtuple(f'{stn}_eec', 'tn_eec', 'ptjet:ptjetch:ptlead:ptleadch:dr:eec')
    self.tn_eec_ch = ROOT.TNtuple(f'{stn}_eec_ch', 'tn_eec_ch', 'ptjet:ptjetch:ptlead:ptleadch:dr:eec')
    
  def get_charged_particles_D0(self, event, ignore_Dstar = True):
    D0_list = []
    D0_list_idx = []
    plist = []
    for ip, p in enumerate([_p for _p in event if _p.isFinal() and _p.isVisible() and _p.isCharged()]):
      if not self.ch_hadron_selector(p):
        continue
      # if p is a charged pion or a charged kaon check if it is a daughter of a D0
      if abs(p.id()) == 211 or abs(p.id()) == 321:
        D0 = is_D0_daughter(p, event)
        if D0:
          if ignore_Dstar:
            if is_Dstar_daughter(D0, event):
              continue
          if not self.D0_selector(D0):
            continue
          if D0.index() in D0_list_idx:
            continue
          D0_list_idx.append(D0.index())
          D0_list.append(D0)
          plist.append(D0)
        else:
          plist.append(p)
      else:
        plist.append(p)
    if len(D0_list) == 0:
      return False
    _h = vector[fj.PseudoJet]()
    for p in plist:
      psj = fj.PseudoJet(p.px(), p.py(), p.pz(), p.e())
      if abs(p.id()) == 421:
        psj.set_user_index(2)
      psj.set_user_index(1)
      _h.push_back(psj)
    h_selected = self.ch_hadron_selector(plist)
    if len(h_selected) == 0:
      return False
    idxs = [psj.user_index() for psj in h_selected]
    if 2 not in idxs:
      return False
    return h_selected

  def analyze(self, pythia, rf):
      # get partons in the last stage before hadronization
      partons = vector[fj.PseudoJet]([fj.PseudoJet(p.px(), p.py(), p.pz(), p.e()) for p in pythia.event if (p.isParton() and (p.isFinal() or (71 <= abs(p.status()) <= 79)))])
      # get hadrons
      hadrons = vector[fj.PseudoJet]([fj.PseudoJet(p.px(), p.py(), p.pz(), p.e()) for p in pythia.event if p.isFinal() and p.isVisible()])
      ch_hadrons = vector[fj.PseudoJet]([fj.PseudoJet(p.px(), p.py(), p.pz(), p.e()) for p in pythia.event if p.isFinal() and p.isVisible() and p.isCharged()])
      # get ch hadrons with daughters of the D0 replaced by the D0
      ch_hadrons_D0 = self.get_charged_particles_D0(pythia.event)
      if ch_hadrons_D0 is False:
        return False
      
      # now run the jet finders
      jet_def_parton = fj.JetDefinition(fj.antikt_algorithm, self.config.jet_R0)
      jet_def_hadron = fj.JetDefinition(fj.antikt_algorithm, self.config.jet_R0)
      jet_selector = fj.SelectorPtMin(self.config.jet_pt_min) * fj.SelectorAbsEtaMax(self.config.jet_abs_eta_max)

      jets_hadron_D0 = fj.sorted_by_pt(self.jet_selector(self.jet_def_hadron(ch_hadrons_D0)))
      if len(jets_hadron_D0) == 0:
        return False
      
      D0_jets = []
      for j in jets_hadron_D0:
        c_idxs = [p.user_index() for p in j.constituents()]
        if 2 not in c_idxs:
          continue
        D0_jets.append(j)
        
      if len(D0_jets) == 0:
        return False
      
      # find jets on parton level
      jets_parton = fj.sorted_by_pt(jet_selector(jet_def_parton(partons)))
      
      # fill things for matched jets
      # self.tn_jet = ROOT.TNtuple(f'{stn}_jet', 'tn_jet', 'ptjet:ptjetch:ptlead:ptleadch')
      # self.tn_eec = ROOT.TNtuple(f'{stn}_eec', 'tn_eec', 'ptjet:ptjetch:ptlead:ptleadch:dr:eec')
      # self.tn_eec_ch = ROOT.TNtuple(f'{stn}_eec_ch', 'tn_eec_ch', 'ptjet:ptjetch:ptlead:ptleadch:dr:eec')

      jet_pairs = []
      for j in jets_parton:
        for jD0 in D0_jets:
          dr = j.delta_R(jD0)
          if dr < self.config.jet_R0:
            jet_pairs.append([j, jD0])

      if len(jet_pairs) == 0:
        return False
      
      # fill the ntuples
      for j, jD0 in jet_pairs:
        ptlead = fj.sorted_by_pt(j.constituents())[0].perp()
        ptleadch = fj.sorted_by_pt(jD0.constituents())[0].perp()
        self.tn_jet.Fill(j.perp(), jD0.perp(), ptlead, ptleadch)
        # calculate EEC for parton jet
        dr, eec = self.calc_eec(j)
        for d, e in zip(dr, eec):
          self.tn_eec.Fill(j.perp(), jD0.perp(), ptlead, ptleadch, d, e)
        # calculate EEC for particles in the hadron jet
        dr, eec = self.calc_eec(jD0)
        for d, e in zip(dr, eec):
          self.tn_eec_ch.Fill(j.perp(), jD0.perp(), ptlead, ptleadch, d, e)
          
      return len(jet_pairs)
        
  def calc_eec(self, j, ptcut=1.0):
    _parts_cut = [p for p in j.constituents() if p.perp() >= ptcut]
    _pairs = list(itertools.product(_parts_cut, repeat=2))
    dr = []
    eec = []
    for first, second in _pairs:
      dr.append(first.delta_R(second))
      eec.append(first.perp() * second.perp() / pow(j.perp(), 2.)
    return dr, eec

      
    
def main():
    parser = argparse.ArgumentParser(description='read hepmc and analyze eecs', prog=os.path.basename(__file__))
    pyconf.add_standard_pythia_args(parser)
    parser.add_argument('--config', help='configure from yaml', type=str)
    parser.add_argument('--input', help='input file', type=str, default=None)
    parser.add_argument('--hepmc', help='2 or 3? - coupled with --input', type=int, default=3)
    parser.add_argument('--write-config', help='write config to yaml and quit', type=str, default='')
    parser.add_argument('-o','--output', help='root output filename', default='analysis_results.root', type=str)
    parser.add_argument('--D0mode', help='set D0 mode', default=0, type=int)
    parser.add_argument('--vetoDstar', help='veto D*->D0pi', default=False, action='store_true')
    parser.add_argument('--ncounts', help='number of D0-jet counts', default=-1, type=int)
    parser.add_argument('--verbose', help="be verbose", default=False, action='store_true')
    parser.add_argument('--debug', help="write debug things", default=False, action='store_true')
    parser.add_argument('--log', help="write stdout to a log file - output file will be same as -o but .log instead of .root", default=False, action='store_true')
    parser.add_argument('--jet-pt-min', help="minimum pT jet to accept", default=20., type=float)
    parser.add_argument('--jet-pt-max', help="max pT jet to accept", default=-1, type=float)
    parser.add_argument('--D0-pt-min', help="minimum pT D0", default=3., type=float)
    parser.add_argument('--D0-pt-max', help="max pT D0", default=100, type=float)
    parser.add_argument('--jet-abs-eta-max', help="max eta of a jet to accept", default=0.5, type=float)
    parser.add_argument('--part-abs-eta-max', help="max eta of a particle to accept", default=1.0, type=float)
    parser.add_argument('--jet-R', help="jet R", default=0.4, type=float)
    parser.add_argument('--parts-select-charged', help="only charged particles", default=False, action='store_true')
    parser.add_argument('--parton-hadron', help="parton then hadron phase", default=False, action='store_true')
    parser.add_argument('--enable-eec', help="enable EEC calculation", default=False, action='store_true')
    parser.add_argument('--', dest='remainder', nargs=argparse.REMAINDER)

    args = parser.parse_args()

    if args.remainder:
        log.info(f'[i] remainder: {args.remainder}')

    config = ConfigData(args=args)
    if args.write_config:
        config.write_to_yaml(args.write_config)
        return

    if args.debug:
        log.set_level(logging.DEBUG)

    if args.output == 'eec_pythia8.root':
        if args.py_vincia:
            args.output = args.output.replace('.root', '_vincia.root')
        if args.py_dire:
            args.output = args.output.replace('.root', '_dire.root')
        log.info(f"[w] using [modified] default output file: {args.output}")
    else:
        log.info(f"[w] using specified output file: {args.output}")

    # abs_eta_max_D0 = 0.8
    abs_rap_max_D0 = 0.8
    # D0_selector = fj.SelectorAbsEtaMax(abs_eta_max_D0) * fj.SelectorPtMin(config.D0_pt_min) * fj.SelectorPtMax(config.D0_pt_max)
    D0_selector = fj.SelectorAbsRapMax(abs_rap_max_D0) * fj.SelectorPtMin(config.D0_pt_min) * fj.SelectorPtMax(config.D0_pt_max)
    if config.D0mode:
        log.info(f'D0 selector: {D0_selector.description()}')

    source_type = 'hepmc'
    pythia = None
    hepmc_input = config.input
    hepmc_event = None 
    
    mycfg = []
    pythia = pyconf.create_and_init_pythia_from_args(config, mycfg)
    if not pythia:
      log.error("[e] pythia initialization failed.")
      return

    if config.nev <= 0 and config.ncounts <= 0:
        config.nev = 10
        log.info(f"[w] setting nev to {config.nev}")
    log.info(f"config: {config}")
    
    fj.ClusterSequence.print_banner()
    print()
    jet_R0 = 0.4
    jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
    max_abs_eta_jet = 0.5
    # jet_selector_ch = fj.SelectorPtMin(args.jet_pt_min) * fj.SelectorPtMax(args.jet_pt_max) * fj.SelectorAbsEtaMax(abs_eta_max_D0 - jet_R0 * 1.02)
    # jet_selector = fj.SelectorPtMin(args.jet_pt_min) * fj.SelectorAbsEtaMax(abs_eta_max_D0 - jet_R0 * 1.02)
    jet_selector_ch = fj.SelectorPtMin(args.jet_pt_min) * fj.SelectorPtMax(args.jet_pt_max) * fj.SelectorAbsEtaMax(max_abs_eta_jet)
    jet_selector = fj.SelectorPtMin(args.jet_pt_min)

    #abs_rap_max_part = 0.9
    #part_selector = fj.SelectorAbsRapMax(abs_rap_max_part) * fj.SelectorPtMin(0.15) * fj.SelectorPtMax(1000.)
    abs_eta_max_part = 0.9
    part_selector = fj.SelectorAbsEtaMax(abs_eta_max_part) * fj.SelectorPtMin(0.15) * fj.SelectorPtMax(1000.)
    log.info(f'part selector: {part_selector.description()}') 

    ptcut = 1.0
    rf = SingleRootFile(fname=config.output)
    rf.root_file.cd()
    tn_jet = ROOT.TNtuple('tn_jet', 'tn_jet', 'ptjet:ptjetch:ptlead:ptleadch')
    tn_eec = ROOT.TNtuple('tn_eec', 'tn_eec', 'ptjet:ptjetch:ptlead:ptleadch:dr:eec')
    tn_eec_ch = ROOT.TNtuple('tn_eec_ch', 'tn_eec_ch', 'ptjet:ptjetch:ptlead:ptleadch:dr:eec')

    _stop = False
    pbar = tqdm.tqdm(range(config.nev))
    pbar.set_description('[i] events processed')
    pbar_counts = tqdm.tqdm(range(config.ncounts))
    pbar_counts.set_description('[i]  counts accepted')
    n_count = 0
    data_source = None
    while not _stop:
        if pythia:
            if not pythia.next():
                continue
            pbar.update(1)
            data_source = pythia
        rv = False
        # fjparts = vector[fj.PseudoJet]([fj.PseudoJet(p.px(), p.py(), p.pz(), p.e()) for p in pythia.event if p.isFinal() and p.isCharged()])
        # fjparts = vector[fj.PseudoJet]([fj.PseudoJet(p.px(), p.py(), p.pz(), p.e()) for p in pythia.event if p.isFinal()])

        fjparts = vector[fj.PseudoJet]()
        log.debug(f'NEW EVENT')

        idx_to_skip = []

        D0list = []    
        if config.D0mode >= 1:
          for ip, p in enumerate(pythia.event):
            # log.debug(f' check: {ip}=?={p.index()}')
            # if p is a D0 and we are in D0 mode
            if abs(p.id()) == 421:
              # if p is a daughter of D* skip it
              Dstar = is_Dstar_daughter(p, pythia.event)
              if Dstar and config.vetoDstar:
                log.debug(f'->-D0: {p.id()} but a D* daughter')
                log.debug(f'     : id={p.id()} pT={p.pT():.3f} eta={p.eta():.3f} y={p.y():.3f} phi={p.phi():.3f} m={p.m():.3f}')
                log.debug(f'   D*: id={Dstar.id()} pT={Dstar.pT():.3f} eta={Dstar.eta():.3f} y={Dstar.y():.3f} phi={Dstar.phi():.3f} m={Dstar.m():.3f}')
                continue
              if config.D0mode == 3 or config.D0mode == 4:
                Kpi = is_D0_Kpi_decay(p, pythia.event)
                if Kpi:
                  log.debug(f'->+D0: {p.id()} a Kpi decay')
                  log.debug(f'     : id={p.id()} pT={p.pT():.3f} eta={p.eta():.3f} y={p.y():.3f} phi={p.phi():.3f} m={p.m():.3f}')
                  log.debug(f'     K: id={pythia.event[Kpi[0]].id()} pT={pythia.event[Kpi[0]].pT():.3f} eta={pythia.event[Kpi[0]].eta():.3f} y={pythia.event[Kpi[0]].y():.3f} phi={pythia.event[Kpi[0]].phi():.3f} m={pythia.event[Kpi[0]].m():.3f}')
                  log.debug(f'     pi: id={pythia.event[Kpi[1]].id()} pT={pythia.event[Kpi[1]].pT():.3f} eta={pythia.event[Kpi[1]].eta():.3f} y={pythia.event[Kpi[1]].y():.3f} phi={pythia.event[Kpi[1]].phi():.3f} m={pythia.event[Kpi[1]].m():.3f}')
                  if config.D0mode == 4:
                    idx_to_skip.extend(Kpi)
                else:
                  log.debug(f'->-D0: {p.id()} but not a Kpi decay')
                  log.debug(f'     : id={p.id()} pT={p.pT():.3f} eta={p.eta():.3f} y={p.y():.3f} phi={p.phi():.3f} m={p.m():.3f}')
                  continue
              psj = fj.PseudoJet(p.px(), p.py(), p.pz(), p.e())
              psj.set_user_index(2)
              # print(f'D0: {psj.perp()} {D0_selector(psj)}')
              if D0_selector(psj):
                if config.D0mode == 2:
                  fjparts.push_back(psj)
                  idx_to_skip.append(ip)
                  _ds = D0_any_decay(p, pythia.event)
                  if _ds:
                    idx_to_skip.extend(_ds)
                log.debug(f'+>+D0: id={p.id()} uidx={psj.user_index()} pT={p.pT():.3f} eta={p.eta():.3f} y={p.y():.3f} phi={p.phi():.3f} m={p.m():.3f}')
                D0list.append(p)
          if len(D0list) == 0:
            continue

        log.debug(f'n of D0s accepted: {len(D0list)}')
        # check the level of the log and if debug print the D0s
        if log.level >= logging.DEBUG:
          for d0 in D0list:
            log.debug(f'  D0: id={d0.id()} pT={d0.pT():.3f} eta={d0.eta():.3f} y={d0.y():.3f} phi={d0.phi():.3f} m={d0.m():.3f}')

        log.debug(f'n parts in vector: {fjparts.size()}')
            
        for ip, p in enumerate(pythia.event):
          if ip in idx_to_skip:
            log.debug(f' -part: id={p.id()} in a list to skip')
            log.debug(f'      : id={p.id()} pT={p.pT():.3f} eta={p.eta():.3f} y={p.y():.3f} phi={p.phi():.3f} m={p.m():.3f}')
            continue
          # if p.isFinal() and p.isVisible():
          # allow neutrinos
          if p.isFinal():
            # if p is a daughter of a D0 skip it
            if config.D0mode == 2 or config.D0mode == 4: 
              D0d = is_D0_daughter(p, pythia.event)
              if D0d:
                log.debug(f' >-part: id={p.id()} but a D0 daughter and mode is {config.D0mode}') 
                log.debug(f'       : id={p.id()} uidx={psj.user_index()} pT={p.pT():.3f} eta={p.eta():.3f} y={p.y():.3f} phi={p.phi():.3f} m={p.m():.3f}')
                log.debug(f'     D0: id={D0d.id()} pT={D0d.pT():.3f} eta={D0d.eta():.3f} y={D0d.y():.3f} phi={D0d.phi():.3f} m={D0d.m():.3f}')
                continue
              # # we do not need this
              # if is_Dstar_daughter(p, pythia.event):
              #   log.debug(f'  -part: {p.id()} but a D* daughter and mode is {config.D0mode}') 
              #   continue
            psj = fj.PseudoJet(p.px(), p.py(), p.pz(), p.e())
            if p.isCharged():
              psj.set_user_index(1)
            else:
              psj.set_user_index(0)
            if part_selector(psj):
              fjparts.push_back(psj)
              log.debug(f'  +part: id={p.id()} uidx={psj.user_index()} pT={p.pT():.3f} eta={p.eta():.3f} y={p.y():.3f} phi={p.phi():.3f} m={p.m():.3f}')
            else:
              # log.debug(f'  -part: id={p.id()} uidx={psj.user_index()} pT={p.pT():.3f} eta={p.eta():.3f} y={p.y():.3f} phi={p.phi():.3f} m={p.m():.3f}')
              pass

        jets = fj.sorted_by_pt(jet_selector(jet_def(fjparts)))
        if len(jets) == 0:
          continue

        _info = Pythia8.getInfo(pythia)
        sigmaGen = _info.sigmaGen()
        ev_weight = _info.weight()
        accepted = False
        for j in jets:
          ptlead = fj.sorted_by_pt(j.constituents())[0].perp()
          _parts_cut = [p for p in j.constituents() if p.perp() >= ptcut]
          _pairs = list(itertools.product(_parts_cut, repeat=2))

          # charged particle jet
          _parts_ch = vector[fj.PseudoJet]()
          _ = [_parts_ch.push_back(p) for p in j.constituents() if p.user_index() > 0]
          
          jetchperp = 0
          ptleadch = 0
          if len(_parts_ch) > 0:
            jetsch = jet_selector_ch(jet_def(_parts_ch))
            if len(jetsch) > 0:
              jetch = jetsch[0]
            else:
              continue
            jetchperp = jetch.perp()
            ptleadch = fj.sorted_by_pt(jetch.constituents())[0].perp()
          else:
            continue

          rv = True
          tn_jet.Fill(j.perp(), jetchperp, ptlead, ptleadch)

          # fill the jet pairs
          for first, second in _pairs:
              dr = first.delta_R(second)
              eec = first.perp() * second.perp() / pow(j.perp(), 2.)
              tn_eec.Fill(j.perp(), jetchperp, ptlead, ptleadch, dr, eec)

          # fill the ch jet pairs
          # _parts_ch_cut = [p for p in j.constituents() if p.user_index() == 0 and p.perp() >= ptcut]
          # _parts_ch_cut = [p for p in j.constituents() if p.user_index() == 1 and p.perp() >= ptcut]
          # _parts_ch_cut = [p for p in j.constituents() if p.user_index() > 0 and p.perp() >= ptcut]
          _parts_ch_cut = [p for p in jetch.constituents() if p.perp() >= ptcut]
          _pairs_ch = list(itertools.product(_parts_ch_cut, repeat=2))
          for first, second in _pairs_ch:
              dr = first.delta_R(second)
              eec = first.perp() * second.perp() / pow(jetchperp, 2.)
              tn_eec_ch.Fill(j.perp(), jetchperp, ptlead, ptleadch, dr, eec)
              
          log.debug(f'* jet: {j.perp()} {j.eta()} {j.phi()}')
          log.debug(f'  n of constituents: {len(j.constituents())} n of cut parts: {len(_parts_cut)} => n of pairs: {len(_pairs)}')
          log.debug(f' +jet: {jetch.perp()} {jetch.eta()} {jetch.phi()}')
          log.debug(f'  n of constituents: {len(jetch.constituents())} n of cut parts: {len(_parts_ch_cut)} => n of pairs: {len(_pairs_ch)}')

        if rv:
            pbar_counts.update(1)
        if pbar.n >= config.nev and config.nev > 0:
            _stop = True
        if pbar_counts.n >= config.ncounts and config.ncounts > 0:
            _stop = True
            
    pbar.close()
    pbar_counts.close()
    if pythia:
        pythia.stat()
        del pythia
    
    rf.close()

if __name__ == '__main__':
    main()