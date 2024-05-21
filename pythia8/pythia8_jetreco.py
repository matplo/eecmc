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

from heppyy.pythia_util import configuration as pyconf
from heppyy.util.logger import Logger
log = Logger()

sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, 'analysis'))
from analysis import ConfigData, JetAnalysis, EECAnalysis, JetChargedFullAnalysis

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


def main():
    parser = argparse.ArgumentParser(description='read hepmc and analyze eecs', prog=os.path.basename(__file__))
    pyconf.add_standard_pythia_args(parser)
    parser.add_argument('--config', help='configure from yaml', type=str)
    parser.add_argument('--input', help='input file', type=str, default=None)
    parser.add_argument('--write-config', help='write config to yaml and quit', type=str, default='')
    parser.add_argument('-o','--output', help='root output filename', default='analysis_results.root', type=str)
    parser.add_argument('--D0mode', help='set D0 mode', default=0, type=int)
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
    if config.nev <= 0 and config.ncounts <= 0:
        config.nev = 10
        log.info(f"[w] setting nev to {config.nev}")
    log.info(f"config: {config}")

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

    D0_selector = fj.SelectorAbsEtaMax(0.8) * fj.SelectorPtMin(config.D0_pt_min) * fj.SelectorPtMax(config.D0_pt_max)
    if config.D0mode:
        log.info(f'[i] D0 selector: {D0_selector.description()}')

    pythia = None
    hepmc = args.input
    if config.input is None:
        mycfg = []
        pythia = pyconf.create_and_init_pythia_from_args(config, mycfg)
        if not pythia:
            print("[e] pythia initialization failed.")
            return
    else:
        print("[w] reading from file not implemented yet")
        return

    # jet_all = JetAnalysis(config, name='jet_all')
    analyses = []

    # jet_parton = None
    # if config.parton_hadron:
    #     jet_parton = JetAnalysis(config, name='jet_parton')
    #     jet_parton.select_parton(True)
    #     analyses.append(jet_parton)

    jfch = JetChargedFullAnalysis(config, name='jfch')
    analyses.append(jfch)

    eec_all = None
    eec_ch = None
    if config.enable_eec:
        eec_all = EECAnalysis(config, name='eec_all', jet_analysis=jfch.j_full_ana)
        analyses.append(eec_all)
        eec_ch = EECAnalysis(config, name='eec_ch', jet_analysis=jfch.j_ch_ana)
        analyses.append(eec_ch)

    _stop = False
    pbar = tqdm.tqdm(range(config.nev))
    pbar.set_description('[i] events processed')
    pbar_counts = tqdm.tqdm(range(config.ncounts))
    pbar_counts.set_description('[i]  counts accepted')
    n_count = 0
    hepmc_event = None
    while not _stop:
        if pythia:
            if not pythia.next():
                continue
        if hepmc:
            hepmc_event = None
        rv = True
        for an in analyses:
            rv = rv and an.analyze_event(pythia=pythia, hepmc=hepmc_event)
        if not rv:
            continue
        pbar.update(1)
        pbar_counts.update(1)
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
