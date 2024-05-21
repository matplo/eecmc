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
from analysis import ConfigData, JetAnalysis, EECAnalysis, JetChargedFullAnalysis, DataSource, SingleRootFile
import hepmc_count_events

def logbins(xmin, xmax, nbins):
        lspace = np.logspace(np.log10(xmin), np.log10(xmax), nbins+1)
        arr = array.array('f', lspace)
        return arr

def main():
    parser = argparse.ArgumentParser(description='read hepmc and analyze eecs', prog=os.path.basename(__file__))
    pyconf.add_standard_pythia_args(parser)
    parser.add_argument('--config', help='configure from yaml', type=str)
    parser.add_argument('--input', help='input file', type=str, default=None)
    parser.add_argument('--hepmc', help='2 or 3? - coupled with --input', type=int, default=3)
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
        log.info(f"[w] using [modified] default output file: {args.output}")
    else:
        log.info(f"[w] using specified output file: {args.output}")

    D0_selector = fj.SelectorAbsEtaMax(0.8) * fj.SelectorPtMin(config.D0_pt_min) * fj.SelectorPtMax(config.D0_pt_max)
    if config.D0mode:
        log.info(f'[i] D0 selector: {D0_selector.description()}')

    source_type = 'hepmc'
    pythia = None
    hepmc = config.input
    hepmc_event = None 
    if config.input is None:
        mycfg = []
        pythia = pyconf.create_and_init_pythia_from_args(config, mycfg)
        source_type = 'pythia'
        if not pythia:
            log.error("[e] pythia initialization failed.")
            return
    else:
        if config.nev < 0:
            config.nev = hepmc_count_events.get_n_per_file(config.input, config.hepmc)[0]
        ###
        # now lets read the HEPMC file and do some jet finding
        if config.hepmc == 3:
            input_hepmc = pyhepmc.io.ReaderAscii(config.input)
        if config.hepmc == 2:
            input_hepmc = pyhepmc.io.ReaderAsciiHepMC2(config.input)
        if input_hepmc.failed():
            log.error("[error] unable to read from {}".format(config.input))
            sys.exit(1)
            return
        hepmc_event = pyhepmc.GenEvent()


    # jet_all = JetAnalysis(config, name='jet_all')
    analyses = []

    # jet_parton = None
    # if config.parton_hadron:
    #     jet_parton = JetAnalysis(config, name='jet_parton')
    #     jet_parton.select_parton(True)
    #     analyses.append(jet_parton)

    rf = SingleRootFile(fname=config.output)
    
    data = DataSource(config, name=f'{source_type}', source_type=source_type)
    data.select_final(True)
    analyses.append(data)

    data_ch = DataSource(config, name=f'{source_type}_ch', source_type=source_type)
    data_ch.select_final(True)
    data_ch.select_charged(True)
    analyses.append(data_ch)

    jet_an = JetAnalysis(config, name='jets', data_source=data)
    analyses.append(jet_an)

    jet_an_ch = JetAnalysis(config, name='jets_ch', data_source=data_ch)
    analyses.append(jet_an_ch)

    eec_all = None
    eec_ch = None
    if config.enable_eec:
        eec_all = EECAnalysis(config, name='eec_all', jet_analysis=jet_an, data_source=data)
        analyses.append(eec_all)
        eec_ch = EECAnalysis(config, name='eec_ch', jet_analysis=jet_an_ch, data_source=data_ch)
        analyses.append(eec_ch)

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
            data_source = pythia
        if hepmc:
            if not input_hepmc.failed():
                _ = input_hepmc.read_event(hepmc_event)
                data_source = hepmc_event
            else:
                _stop = True
                break
        rv = True
        for an in analyses:
            rv = rv and an.process_event(data_source)
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
    if pythia:
        pythia.stat()
        del pythia
    
    rf.close()

if __name__ == '__main__':
    main()
