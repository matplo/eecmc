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
from analysis import ConfigData, SingleRootFile

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

                
def calc_eec(j, ptcut=1.0):
  _parts_cut = [p for p in j.constituents() if p.perp() >= ptcut]
  _pairs = list(itertools.product(_parts_cut, repeat=2))
  dr = []
  eec = []
  for first, second in _pairs:
    dr.append(first.delta_R(second))
    eec.append(first.perp() * second.perp() / pow(j.perp(), 2.))
  return dr, eec

def calc_eec_list(l, ptcut=1.0, weight=1.0):
  _parts_cut = [p for p in l if p.perp() >= ptcut]
  _pairs = list(itertools.product(_parts_cut, repeat=2))
  dr = []
  eec = []
  for first, second in _pairs:
    dr.append(first.delta_R(second))
    eec.append(first.perp() * second.perp() / pow(weight, 2.))
  return dr, eec

# D0select:
# 0 - no D0 selection
# 1 - D0 selection - but not a D* daughter
# 2 - D0 selection and remove D0 daughters
# 4 - D0 selection Kpi decay and remove D0 daughters

anLevelTypes = {
    'parton': 0,
    'hadron': 1,
    'charged': 2,
    'D0': 3,
    'D0Kpi': 4,
}

def get_parts_from_pythia(pythia):
  d_fjparts = {}
  for iltype, ltype in enumerate(anLevelTypes):
    fjparts = vector[fj.PseudoJet]()
    D0_indexes_ltype_3 = []
    D0_indexes_ltype_4 = []
    for p in pythia.event:
      if iltype == 0:
        if (p.isParton()) and ( (p.isFinal()) or (71 <= abs(p.status()) and abs(p.status()) <= 79) ):
          psj = fj.PseudoJet(p.px(), p.py(), p.pz(), p.e())
          psj.set_user_index(p.index())      
          fjparts.push_back(psj)
      if p.isParton():
        continue
      if iltype == 1:
        if p.isFinal() and p.isVisible():
          psj = fj.PseudoJet(p.px(), p.py(), p.pz(), p.e())
          psj.set_user_index(p.index())      
          fjparts.push_back(psj)
      if iltype == 2:
        if p.isFinal() and p.isVisible() and p.isCharged():
          psj = fj.PseudoJet(p.px(), p.py(), p.pz(), p.e())
          psj.set_user_index(p.index())      
          fjparts.push_back(psj)
      if iltype == 3:
        if p.isFinal() and p.isVisible() and p.isCharged():
          D0 = is_D0_daughter(p, pythia.event)
          Dstar = is_Dstar_daughter(p, pythia.event)
          if Dstar:
            continue
          if D0:
            if D0.index() in D0_indexes_ltype_3:
              continue
            psj = fj.PseudoJet(D0.px(), D0.py(), D0.pz(), D0.e())
            D0_indexes_ltype_3.append(D0.index())
            psj.set_user_index(p.index())      
            fjparts.push_back(psj)
          else:
            psj = fj.PseudoJet(p.px(), p.py(), p.pz(), p.e())
            psj.set_user_index(p.index())      
            fjparts.push_back(psj)
      if iltype == 4:
        if p.isFinal() and p.isVisible() and p.isCharged():
          D0 = is_D0_daughter(p, pythia.event)
          Dstar = is_Dstar_daughter(p, pythia.event)
          if Dstar:
            continue
          if D0 and (abs(p.id()) == 321 or abs(p.id()) == 211):
            if D0.index() in D0_indexes_ltype_4:
              continue
            psj = fj.PseudoJet(D0.px(), D0.py(), D0.pz(), D0.e())
            D0_indexes_ltype_4.append(D0.index())
            psj.set_user_index(p.index())      
            fjparts.push_back(psj)
          else:
            psj = fj.PseudoJet(p.px(), p.py(), p.pz(), p.e())
            psj.set_user_index(p.index())      
            fjparts.push_back(psj)
    if iltype < 3:
      d_fjparts[ltype] = fjparts
    else:
      if iltype == 3 and len(D0_indexes_ltype_3) > 0:
        d_fjparts[ltype] = fjparts
      if iltype == 4 and len(D0_indexes_ltype_4) > 0:
        d_fjparts[ltype] = fjparts
  return d_fjparts

################################################################################

class EECAnalysis(object):
  def __init__(self, ltype, config):
    self.ltype = ltype
    self.tn_jet = ROOT.TNtuple(f'tn_jet_{ltype}', 'tn_jet', 'ptjet:ptjetch:ptlead:ptleadch:etajet')
    self.tn_eec = ROOT.TNtuple(f'tn_jet_{ltype}_eec', 'tn_eec', 'ptjet:ptjetch:ptlead:ptleadch:dr:eec:etajet')
    self.tn_eec_ch = ROOT.TNtuple(f'tn_jet_ch_{ltype}_eec', 'tn_eec_ch', 'ptjet:ptjetch:ptlead:ptleadch:dr:eec:etajet')

  def fill_jets(self, jets, pythia):
    for j in jets:
      ptlead = fj.sorted_by_pt(j.constituents())[0].perp()
      _ch_parts = [p for p in j.constituents() if pythia.event[p.user_index()].isCharged()]
      ptjetch = sum([p.perp() for p in _ch_parts])
      if len(_ch_parts) > 0:
        ptleadch = fj.sorted_by_pt(_ch_parts)[0].perp()
      else:
        ptleadch = 0
      self.tn_jet.Fill(j.perp(), ptjetch, ptlead, ptleadch, j.eta())
      dr, eec = calc_eec(j, ptcut=1.0)
      for d, e in zip(dr, eec):
        self.tn_eec.Fill(j.perp(), ptjetch, ptlead, ptleadch, d, e, j.eta())
      dr, eec = calc_eec_list(_ch_parts, ptcut=1.0, weight=j.perp())
      for d, e in zip(dr, eec):
        self.tn_eec_ch.Fill(j.perp(), ptjetch, ptlead, ptleadch, d, e, j.eta())
    
################################################################################

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

    # ARGUMENTS ################################################################################
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

    # SELECTORS ################################################################################

    ptcut = 1.0

    # abs_eta_max_D0 = 0.8
    abs_rap_max_D0 = 0.8
    # D0_selector = fj.SelectorAbsEtaMax(abs_eta_max_D0) * fj.SelectorPtMin(config.D0_pt_min) * fj.SelectorPtMax(config.D0_pt_max)
    D0_selector = fj.SelectorAbsRapMax(abs_rap_max_D0) * fj.SelectorPtMin(config.D0_pt_min) * fj.SelectorPtMax(config.D0_pt_max)
    if config.D0mode:
        log.info(f'D0 selector: {D0_selector.description()}')

    fj.ClusterSequence.print_banner()
    print()
    jet_R0 = 0.4
    jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
    max_abs_eta_jet = 0.5
    jet_selector_ch = fj.SelectorPtMin(config.jet_pt_min) * fj.SelectorPtMax(config.jet_pt_max) * fj.SelectorAbsEtaMax(max_abs_eta_jet)
    jet_selector = fj.SelectorPtMin(config.jet_pt_min)
    jet_selector_D0 = fj.SelectorPtMin(config.jet_pt_min) * fj.SelectorPtMax(config.jet_pt_max) * fj.SelectorAbsEtaMax(max_abs_eta_jet)

    #abs_rap_max_part = 0.9
    #part_selector = fj.SelectorAbsRapMax(abs_rap_max_part) * fj.SelectorPtMin(0.15) * fj.SelectorPtMax(1000.)
    abs_eta_max_part = 0.9
    part_selector = fj.SelectorAbsEtaMax(abs_eta_max_part) * fj.SelectorPtMin(0.15) * fj.SelectorPtMax(1000.)
    log.info(f'part selector: {part_selector.description()}') 

    # PYTHIA CONFIG ################################################################################
    
    source_type = 'hepmc'
    pythia = None
    hepmc_input = config.input
    hepmc_event = None 

    mycfg = []
    pythia = pyconf.create_and_init_pythia_from_args(config, mycfg)
    if not pythia:
      log.error("[e] pythia initialization failed.")
      return

    # EVENT CONFIG ################################################################################

    if config.nev <= 0 and config.ncounts <= 0:
        config.nev = 10
        log.info(f"[w] setting nev to {config.nev}")
    log.info(f"config: {config}")
    
    # OUTPUT FILE CONFIG ################################################################################
    rf = SingleRootFile(fname=config.output)
    rf.root_file.cd()
    d_jan = {}
    for ltype in anLevelTypes:
      d_jan[ltype] = EECAnalysis(ltype, config)

    # EVENT LOOP ################################################################################
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
      d_parts = get_parts_from_pythia(pythia)
      d_jets = {}
      for iltype, ltype in enumerate(d_parts):
        log.debug(f'type={iltype} - {ltype} {len(d_parts[ltype])}')
        jets = fj.sorted_by_pt(jet_selector(jet_def(d_parts[ltype])))
        if len(jets) == 0:
          continue
        log.debug(f'jets: {len(jets)}')
        d_jets[ltype] = jets
      log.debug(f'------')

      if 'D0Kpi' in d_jets:
        if len(d_jets['D0Kpi']) > 0:
          rv = 1

      for ltype in d_jets:
        d_jan[ltype].fill_jets(d_jets[ltype], pythia)
          
      # analyze EECs from the jets
      # write ntuple to a file...    
      
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