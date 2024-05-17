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
            log.warning('overriding values from command line')
            for s in sys.argv:
                if s.startswith('--'):
                    s = s[2:].replace('-', '_')
                    for k, v in self.args.__dict__.items():
                        if s != k:
                            continue
                        log.warning(f'overriding {k} with {v}')
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
        self._log = Logger()
        self.config_input = config
        self.cfg = ConfigData()
        if isinstance(config, str):
            self.cfg.configure_from_yaml(config)
        if isinstance(config, dict):
            self.cfg.configure_from_dict(config)
        if isinstance(config, ConfigData):
            self.cfg.configure_from_dict(config.__dict__)
        super(AnalysisBase, self).__init__(**kwargs)
        if self.name is None:
            self.name = 'AnalysisBase'
        self.name = NamesStore().add_unique(self.name)
        self._nev = 0
        self.open_root_file()
        self.init()
        self._log.info(f"{self.name}: initialized - with config {self.cfg.__dict__}")

    def open_root_file(self):
        # open a root file if needed
        if self.root_file is None:
            self.objects = []
            if self.cfg.output is None:
                self.cfg.output = 'analysis_output.root'
                self.cfg.output.replace('.root', f'_{self.name}.root')
            if self.cfg.output.endswith('.root'):
                self.cfg.output = self.cfg.output.replace('.root', f'_{self.name}.root')
            self.cfg.output = NamesStore().add_unique(self.cfg.output)
            if self.root_file is None:
                self.root_file = ROOT.TFile(self.cfg.output, 'RECREATE')
                self._log.info(f'{self.name}: opened root file {self.root_file.GetName()}')
        self.root_file.cd()

    def init(self):
        self._log.debug('AnalysisBase::init()')

    def analyze_event(self, pythia=None, hepmc=None, nev=None):
        if nev is not None:
            self._nev = nev
        else:
            self._nev += 1
        rv = True
        if pythia:
            rv = rv and self.analyze_pythia_event(pythia)
        if hepmc:
            rv = rv and self.analyze_hepmc_event(hepmc)
        return rv

    def analyze_pythia_event(self, pythia):
        pass

    def analyze_hepmc_event(self, hepmc):
        pass

    def __del__(self):
        if self.root_file:
            self.root_file.Write()
            self.root_file.Close()
            self._log.info(f'{self.name}: wrote {self.root_file.GetName()}')
        else:
            self._log.debug('no root file to write')
        self.root_file = None

class JetAnalysis(AnalysisBase):

    def init(self):
        super(JetAnalysis, self).init()
        self.open_root_file()
        self.root_file.cd() # done already in the base class but just to be sure
        self.tn_events 	= ROOT.TNtuple(f'tn_events_{self.name}', 	'tn_events', 	'nev:xsec:ev_weight:njets:nparts')
        self.tn_jet 	= ROOT.TNtuple(f'tn_jet_{self.name}', 		'tn_jet', 		'nev:ij:pt:eta:phi:m:xsec:ev_weight')
        self.tn_parts 	= ROOT.TNtuple(f'tn_parts_{self.name}', 	'tn_parts', 	'nev:idx:ij:pt:eta:phi:m:pid:status:xsec:ev_weight')
        # jet finder
        self.jet_def 			= fj.JetDefinition(fj.antikt_algorithm, self.cfg.jet_R)
        self.jet_selector = fj.SelectorPtMin(self.cfg.jet_pt_min)
        if self.part_abs_eta_max:
            self.part_selector = fj.SelectorAbsEtaMax(self.cfg.part_abs_eta_max)
        if self.cfg.jet_abs_eta_max > 0:
            if (
                self.cfg.part_abs_eta_max
                and self.cfg.jet_abs_eta_max
                > self.cfg.part_abs_eta_max - self.cfg.jet_R
            ):
                self.cfg.jet_abs_eta_max = self.cfg.part_abs_eta_max - self.cfg.jet_R
                self._log.warning(
                    f'{self.name}: jet_abs_eta_max adjusted to {self.cfg.jet_abs_eta_max}'
                )
            self.jet_selector = self.jet_selector * fj.SelectorAbsEtaMax(
                self.cfg.jet_abs_eta_max
            )
        if self.cfg.jet_eta_max is not None:
            self.jet_selector = self.jet_selector * fj.SelectorEtaMax(
                self.cfg.jet_eta_max
            )
        if self.cfg.jet_eta_min is not None:
            self.jet_selector = self.jet_selector * fj.SelectorEtaMin(
                self.cfg.jet_eta_min
            )
        if self.cfg.jet_pt_max > 0:
            self.jet_selector = self.jet_selector * fj.SelectorPtMax(
                self.cfg.jet_pt_max
            )
        self._log.info(f'{self.name}: jet definition: {self.jet_def.description()}')
        self._log.info(f'{self.name}: jet selector: {self.jet_selector.description()}')
        self._log.info(
            f'{self.name}: particle selection set: self.parts_select_final={self.cfg.parts_select_final} self.parts_select_charged={self.cfg.parts_select_charged}'
        )
        if self.cfg.parts_select_parton is None:
            self.select_parton(False)

    def select_final(self, select=True):
        self.cfg.parts_select_final = select
        self._log.info(
            f'{self.name}: particle selection set: self.parts_select_final={self.cfg.parts_select_final} self.parts_select_charged={self.cfg.parts_select_charged}'
        )

    def select_charged(self, select=True):
        self.cfg.parts_select_charged = select
        self._log.info(
            f'{self.name}: particle selection set: self.parts_select_final={self.cfg.parts_select_final} self.parts_select_charged={self.cfg.parts_select_charged}'
        )

    def select_parton(self, select=True):
        self.cfg.parts_select_parton = select
        self._log.info(
            f'{self.name}: particle selection set: self.parts_select_parton={self.cfg.parts_select_parton}'
        )

    def select_pythia_particles(self, pythia):
        self.py_parts = []
        self.parts = vector[fj.PseudoJet]()
        if self.cfg.parts_select_parton:
            self.py_parts = [(i, p, p.id(), p.status()) for i, p in enumerate(pythia.event) if (p.isParton() and (p.isFinal() or (71 <= abs(p.status()) <= 79)))]
            # note: a PID of 1103 represents a diquark made up of an up quark (1) and a down anti-quark (1), with a spin of 1 (the third digit), and a "j" value of 3 (the last digit).
        else:
            if self.cfg.parts_select_final:
                if self.cfg.parts_select_charged:
                    self.py_parts = [(i, p, p.id(), p.status()) for i, p in enumerate(pythia.event) if p.isFinal() and p.isCharged()]
                else:
                    self.py_parts = [(i, p, p.id(), p.status()) for i, p in enumerate(pythia.event) if p.isFinal()]
            else:
                self.py_parts = [(i, p, p.id(), p.status()) for i, p in enumerate(pythia.event)]
                if self._nev < 1:
                    self._log.warning(f'{self.name}: selecting all particles final and not final state!')
        for xp in self.py_parts:
            _psj = fj.PseudoJet(xp[1].px(), xp[1].py(), xp[1].pz(), xp[1].e())
            _psj.set_user_index(xp[0])
            self.parts.push_back(_psj)
        if self.part_selector:
            self.parts = self.part_selector(self.parts)

    def analyze_pythia_event(self, pythia):
        self.select_pythia_particles(pythia)
        if len(self.parts) < 1:
            self._log.debug(f'{self.name}: no particles to analyze len(self.parts)={len(self.parts)}')
            return False
        self.pythia_info = Pythia8.getInfo(pythia)
        self.ev_weight = self.pythia_info.weight()
        self.xsec = self.pythia_info.sigmaGen()
        return self.run_jet_analysis()

    def run_jet_analysis(self):
        self.jets = fj.sorted_by_pt(self.jet_selector(self.jet_def(self.parts)))
        if len(self.jets) < 1:
            return False
        self.tn_events.Fill(self._nev, self.xsec, self.ev_weight, len(self.jets), len(self.parts))
        for ij, j in enumerate(self.jets):
            self.tn_jet.Fill(self._nev, ij, j.pt(), j.eta(), j.phi(), j.m(), self.xsec, self.ev_weight)
        if self.cfg.write_parts:
            for p, py_p in zip(self.parts, self.py_parts):
                ijet = [i for i, j in enumerate(self.jets) if p.user_index() in [x.user_index() for x in j.constituents()]]
                if len(ijet) < 1:
                    ijet = -1
                else:
                    ijet = ijet[0]
                self.tn_parts.Fill(self._nev, p.user_index(), ijet, p.pt(), p.eta(), p.phi(), p.m(), py_p[2], py_p[3], self.xsec, self.ev_weight)
        return True


class JetChargedFullAnalysis(AnalysisBase):

    def init(self):
        super(JetChargedFullAnalysis, self).init()
        self.open_root_file()
        self.root_file.cd()  # done already in the base class but just to be sure
        self.tn_events = ROOT.TNtuple(f"tn_events_{self.name}", "tn_events", "nev:xsec:ev_weight:njets:njets_ch:nparts:nparts_ch")
        self.tn_correl = ROOT.TNtuple(
            f"tn_correl_{self.name}",
            "tn_correl",
            "nev:ij:pt:eta:phi:pt_ch_in_full:xsec:ev_weight:ij_ch:pt_ch:eta_ch:phi_ch:pt_ratio:pt_ratio_ch:n_matched",
        )
        self.tn_correl_parton_full = ROOT.TNtuple(
            f"tn_correl_parton_full_{self.name}",
            "tn_correl_parton_full",
            "nev:ij:pt:eta:phi:pt_h2p:xsec:ev_weight:ij_h:pt_h:eta_h:phi_h:pt_ratio:n_matched",
        )
        self.tn_correl_parton_ch = ROOT.TNtuple(
            f"tn_correl_parton_ch_{self.name}",
            "tn_correl_parton_ch",
            "nev:ij:pt:eta:phi:pt_h2p:xsec:ev_weight:ij_h:pt_h:eta_h:phi_h:pt_ratio:n_matched",
        )
        # note the pt cut on jets will be applied to charged jets - the it is somewhat effectively lower for full jets
        self.j_ch_ana = JetAnalysis(self.cfg, name='jet_charged')
        self.j_ch_ana.select_charged(select=True)
        self.j_ch_ana.select_final(select=True)
        self.j_full_ana = JetAnalysis(self.cfg, name='jet_full')
        self.j_full_ana.select_charged(select=False)
        self.j_full_ana.select_final(select=True)
        if self.cfg.parton_hadron:
            self.j_parton_ana = JetAnalysis(self.cfg, name="jet_parton")
            self.j_parton_ana.select_parton(select=True)

    def analyze_pythia_event(self, pythia):
        rv = True
        if self.cfg.parton_hadron:
            rv = rv and self.j_parton_ana.analyze_event(pythia=pythia, nev=self._nev)
        if not rv:
            return False
        rv = rv and self.j_ch_ana.analyze_event(pythia=pythia, nev=self._nev)
        if not rv:
            return False
        rv = rv and self.j_full_ana.analyze_event(pythia=pythia, nev=self._nev)
        if not rv:
            return False
        self.tn_events.Fill(self._nev, self.j_full_ana.xsec, self.j_full_ana.ev_weight,
                            len(self.j_full_ana.jets), len(self.j_ch_ana.jets),
                            len(self.j_full_ana.parts), len(self.j_ch_ana.parts))
        n_matched = {}
        for ij, j in enumerate(self.j_full_ana.jets):
            sum_charged_in_full = sum([int(pythia.event[p.user_index()].isCharged()) * p.perp() for p in j.constituents()])
            for ijch, jch in enumerate(self.j_ch_ana.jets):
                if j.delta_R(jch) > self.cfg.jet_R / 2.0:
                    continue
                n_matched[ij] = n_matched.get(ij, 0) + 1
                pt_ratio = jch.pt() / j.pt()
                pt_ratio_charged = jch.pt() / sum_charged_in_full
                # x-check:
                self._log.debug(f'{self.name}: pt_ratio={pt_ratio}')
                self._log.debug(f'{self.name}: pt_ratio_charged={pt_ratio_charged}')
                if self._log.logger.getEffectiveLevel() == logging.DEBUG:
                    for p in jch.constituents():
                        if p.user_index() not in [x.user_index() for x in j.constituents()]:
                            self._log.error(f'{self.name}: particle {p.user_index()} in ch jet not found in full jet')
                self.tn_correl.Fill(self._nev, ij, j.pt(), j.eta(), j.phi(), sum_charged_in_full,
                                    self.j_full_ana.pythia_info.sigmaGen(), self.j_full_ana.pythia_info.weight(),
                                    ijch, jch.pt(), jch.eta(), jch.phi(),
                                    pt_ratio, pt_ratio_charged, n_matched[ij])
        if self.cfg.parton_hadron:
            n_matched = {}
            for ij, j in enumerate(self.j_parton_ana.jets):
                for ijch, jch in enumerate(self.j_full_ana.jets):
                    if j.delta_R(jch) > self.cfg.jet_R / 2.0:
                        continue
                    n_matched[ij] = n_matched.get(ij, 0) + 1
                    pt_ratio = jch.pt() / j.pt()
                    self.tn_correl_parton_full.Fill(self._nev, ij, j.pt(), j.eta(), j.phi(), sum_charged_in_full,
                                            self.j_full_ana.pythia_info.sigmaGen(), self.j_full_ana.pythia_info.weight(),
                                            ijch, jch.pt(), jch.eta(), jch.phi(),
                                            pt_ratio, n_matched[ij])
            n_matched = {}
            for ij, j in enumerate(self.j_parton_ana.jets):
                for ijch, jch in enumerate(self.j_ch_ana.jets):
                    if j.delta_R(jch) > self.cfg.jet_R / 2.0:
                        continue
                    n_matched[ij] = n_matched.get(ij, 0) + 1
                    pt_ratio = jch.pt() / j.pt()
                    self.tn_correl_parton_ch.Fill(self._nev, ij, j.pt(), j.eta(), j.phi(), sum_charged_in_full,
                                            self.j_full_ana.pythia_info.sigmaGen(), self.j_full_ana.pythia_info.weight(),
                                            ijch, jch.pt(), jch.eta(), jch.phi(),
                                            pt_ratio, n_matched[ij])

        return rv

class EECAnalysis(AnalysisBase):

    def init(self):
        super(EECAnalysis, self).init()
        self.tn_eec = ROOT.TNtuple(f'tn_eec_{self.name}', 'tn_eec', 'nev:ij:dr:pt1:pt2:eec:ptjet:xsec:ev_weight:ptcut')
        self.pt_cuts = [0, 0.15, 1, 2]
        self.own_jet_analysis = False
        if self.jet_analysis is None:
            self.jet_analysis = JetAnalysis(self.cfg, name=f'{self.name}_jet_analysis')
            self.jet_analysis.select_final(True)
            self.jet_analysis.select_charged(True)
            self.own_jet_analysis = True

    def analyze_pythia_event(self, pythia):
        if self.own_jet_analysis:
            rv = self.jet_analysis.analyze_event(pythia=pythia, nev=self._nev)
            if not rv:
                return False
        rvs = [self.analyze_eec(ptcut) for ptcut in self.pt_cuts]
        return any(rvs)

    def analyze_eec(self, ptcut):
        # Generate all pairs from parts, excluding pairs of the same element
        # self.pairs = list(itertools.combinations(parts, 2))
        # Generate all pairs from parts, including pairs of the same element
        for ij, j in enumerate(self.jet_analysis.jets):
            _parts_cut = [p for p in j.constituents() if p.perp() > ptcut]
            _pairs = list(itertools.product(_parts_cut, repeat=2))
            log.debug(f'{self.name}: number of pairs: {len(_pairs)} with ptcut: {ptcut}')
            if len(_pairs) < 1:
                continue
            for first, second in _pairs:
                dr = first.delta_R(second)
                eec = first.perp() * second.perp() / pow(j.perp(), 2.)
                self.tn_eec.Fill(self._nev, ij, dr, first.perp(), second.perp(), eec, j.perp(), self.jet_analysis.xsec, self.jet_analysis.ev_weight, ptcut)
        return True


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
