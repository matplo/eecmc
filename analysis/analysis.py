import yaml
import sys
import ROOT
import logging
import itertools
import pyhepmc
import particle

import heppyy.util.fastjet_cppyy
import heppyy.util.pythia8_cppyy
import heppyy.util.heppyy_cppyy

from cppyy.gbl import fastjet as fj
from cppyy.gbl import Pythia8
from cppyy.gbl.std import vector

from yasp import GenericObject
from heppyy.util.logger import Logger
log = Logger()

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

# singleton class for for analysis output root file
class SingleRootFile(object):

    __instance = None

    def __new__(cls, fname):
        if SingleRootFile.__instance is None:
            SingleRootFile.__instance = object.__new__(cls)
            SingleRootFile.__instance.filename = fname
            SingleRootFile.__instance.root_file = ROOT.TFile(fname, 'RECREATE')            
            SingleRootFile.__instance.objects = []
            log.info(f'new root file: {SingleRootFile.__instance.root_file.GetName()}')
        return SingleRootFile.__instance

    def add(self, obj):
        if obj:
            self.objects.append(obj)

    def __del__(self):
        if SingleRootFile.__instance:
            if SingleRootFile.__instance.root_file:
                log.info(f'closing {SingleRootFile.__instance.root_file.GetName()}')
                SingleRootFile.__instance.close()
        
    def close(self):
        _rfile = self.root_file
        log.info(f'root_file is {_rfile}')
        if _rfile:
            _rfile.cd()
            _ = [o.Write() for o in self.objects]
            _rfile.Write()
            _rfile.Close()
            log.info(f'wrote {_rfile.GetName()}')
            log.info(f'purging {_rfile.GetName()}')
            _rfile = ROOT.TFile(self.filename, 'UPDATE')
            _rfile.Purge()
            _rfile.Write()
            _rfile.Close()
        else:
            log.debug('no root file to write')
        _rfile = None
        log.info(f'done {self.filename}')

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
        if self.source_type is None and self.data_source is None:
            raise ValueError(f'[e] {self.name} source_type nor data_source not specified')
        if self.name is None:
            self.name = 'AnalysisBase'
        self.name = NamesStore().add_unique(self.name)
        self._nev = 0
        self.open_root_file()
        self.init()
        self._log.debug(f"{self.name}: initialized - with config {self.cfg.__dict__}")

    def open_root_file(self):
        # open a root file if needed
        if self.separate_root_file:
            if self.root_file is None:
                self.objects = []
                if self.cfg.output is None:
                    self.cfg.output = f'analysis_output.root'
                    self.cfg.output.replace('.root', f'_{self.name}.root')
                if self.cfg.output.endswith('.root'):
                    self.cfg.output = self.cfg.output.replace('.root', f'_{self.name}.root')
                self.cfg.output = NamesStore().add_unique(self.cfg.output)
                if self.root_file is None:
                    self.root_file = ROOT.TFile(self.cfg.output, 'RECREATE')
                    self._log.info(f'{self.name}: opened root file {self.root_file.GetName()}')
                    self.owns_root_file = True
        else:
            _srf = SingleRootFile(self.cfg.output)
            _srf.add(self)
            self.root_file = _srf.root_file
            self.owns_root_file = False
        self.root_file.cd()

    def init(self):
        self._log.debug('AnalysisBase::init()')

    def process_event(self, input_object, nev=None):
        if nev is not None:
            self._nev = nev
        else:
            self._nev += 1
        return True

    def Write(self):
        self.finalize()

    def finalize(self):
        if self.root_file:
            self._log.info(f'{self.name}: writing to {self.root_file.GetName()}')
            self.root_file.Write()
            if self.owns_root_file:
                self._log.info(f'{self.name}: closing {self.root_file.GetName()}')
                self.root_file.Close()
        else:
            self._log.debug('no root file to write')
        self.root_file = None
        self.done = True

    def __del__(self):
        if self.done:
            return
        self.finalize()

# simplify this class with a Generator class able to read hepmc events, or pythia, or something else...
# should initialize with an event and type of event
# should read all the particles in as specified by the base class
# Generator subclassess: GenPythia, GenHepMC, etc.
# should hold the event and the particles

class DataSource(AnalysisBase):
    __known_sources__ = ['pythia', 'hepmc']

    def init(self):
        super(DataSource, self).init()
        if self.source_type not in self.__known_sources__:
            raise ValueError(f'[e] unknown data source type {self.source_type}')
        self.open_root_file()
        self.root_file.cd() # done already in the base class but just to be sure
        self.tn_events 	= ROOT.TNtuple(f'tn_events_{self.name}', 'tn_events', 'nev:xsec:ev_weight:nparts')
        if self.cfg.write_parts:
            self.tn_parts 	= ROOT.TNtuple(f'tn_parts_{self.name}', 'tn_parts', 'nev:xsec:ev_weight:idx:pt:eta:phi:m:pid:status')
        if self.part_abs_eta_max:
            self.part_selector = fj.SelectorAbsEtaMax(self.cfg.part_abs_eta_max)
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

    def process_event(self, input_object, nev=None):
        rv = super(DataSource, self).process_event(input_object, nev=nev)
        if self.source_type == 'pythia':
            self.source_object = input_object
            rv = rv and self.process_pythia(input_object)
        if self.source_type == 'hepmc':
            self.source_object = input_object
            rv = rv and self.process_hepmc(input_object)
        return rv

    def select_hepmc_particles(self, hepmc):
        self.source_parts = []
        fjparts = []
        for i,p in enumerate(hepmc.particles):
            if p.status == 1 and not p.end_vertex:
                if self.cfg.parts_select_charged:
                    _charge = None
                    try: 
                        _charge = particle.Particle.from_pdgid(p.pid).charge
                    except particle.particle.particle.ParticleNotFound as e:
                        self._log.debug('particle not found:', e)
                        continue
                    if _charge == 0:
                        continue
                self.source_parts.append((i, p, p.pid, p.status))
                psj = fj.PseudoJet(p.momentum.px, p.momentum.py, p.momentum.pz, p.momentum.e)
                psj.set_user_index(i)
                fjparts.append(psj)
        self.parts = vector[fj.PseudoJet](fjparts)
        if self.part_selector:
            self.parts = self.part_selector(self.parts)
        if len(self.parts) < 1:
            return False
            
    def process_hepmc(self, hepmc):
        self.select_hepmc_particles(hepmc)
        if len(self.parts) < 1:
            self._log.debug(f'{self.name}: no particles to analyze len(self.parts)={len(self.parts)}')
            return False
        self.ev_weight = hepmc.weight()
        self.xsec = hepmc.cross_section.xsec()
        return self.fill_ntuples()

    def select_pythia_particles(self, pythia):
        self.source_parts = []
        self.parts = vector[fj.PseudoJet]()
        if self.cfg.parts_select_parton:
            self.source_parts = [(i, p, p.id(), p.status()) for i, p in enumerate(pythia.event) if (p.isParton() and (p.isFinal() or (71 <= abs(p.status()) <= 79)))]
            # note: a PID of 1103 represents a diquark made up of an up quark (1) and a down anti-quark (1), with a spin of 1 (the third digit), and a "j" value of 3 (the last digit).
        else:
            if self.cfg.parts_select_final:
                if self.cfg.parts_select_charged:
                    self.source_parts = [(i, p, p.id(), p.status()) for i, p in enumerate(pythia.event) if p.isFinal() and p.isCharged()]
                else:
                    self.source_parts = [(i, p, p.id(), p.status()) for i, p in enumerate(pythia.event) if p.isFinal()]
            else:
                self.source_parts = [(i, p, p.id(), p.status()) for i, p in enumerate(pythia.event)]
                if self._nev < 1:
                    self._log.warning(f'{self.name}: selecting all particles final and not final state!')
        for xp in self.source_parts:
            _psj = fj.PseudoJet(xp[1].px(), xp[1].py(), xp[1].pz(), xp[1].e())
            _psj.set_user_index(xp[0])
            self.parts.push_back(_psj)
        if self.part_selector:
            self.parts = self.part_selector(self.parts)
        if len(self.parts) < 1:
            return False

    def process_pythia(self, pythia):
        self.select_pythia_particles(pythia)
        if len(self.parts) < 1:
            self._log.debug(f'{self.name}: no particles to analyze len(self.parts)={len(self.parts)}')
            return False
        self.pythia_info = Pythia8.getInfo(pythia)
        self.ev_weight = self.pythia_info.weight()
        self.xsec = self.pythia_info.sigmaGen()
        return self.fill_ntuples()

    def fill_ntuples(self):
        self.tn_events.Fill(self._nev, self.xsec, self.ev_weight, len(self.parts))
        if self.tn_parts:
            for p, source_p in zip(self.parts, self.source_parts):
                # 'nev:xsec:ev_weight:idx:pt:eta:phi:m:pid:status'
                self.tn_parts.Fill(self._nev, self.xsec, self.ev_weight, p.user_index(), p.pt(), p.eta(), p.phi(), p.m(), source_p[2], source_p[3])
        return True

class JetAnalysis(AnalysisBase):

    def init(self):
        super(JetAnalysis, self).init()
        self.open_root_file()
        self.root_file.cd() # done already in the base class but just to be sure
        self.tn_jet = ROOT.TNtuple(f'tn_jet_{self.name}', 'tn_jet', 'nev:xsec:ev_weight:nj:ij:pt:eta:phi:m')
        # every jet in the file gets a unique int - nj
        self.njet_count = 0
        # jet finder
        self.jet_definition = fj.JetDefinition(fj.antikt_algorithm, self.cfg.jet_R)
        self.jet_selector = fj.SelectorPtMin(self.cfg.jet_pt_min)
        if self.cfg.jet_abs_eta_max > 0:
            if (
                self.data_source.part_abs_eta_max
                and self.cfg.jet_abs_eta_max
                > self.data_source.part_abs_eta_max - self.cfg.jet_R
            ):
                self.cfg.jet_abs_eta_max = self.data_source.part_abs_eta_max - self.cfg.jet_R
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
        self._log.info(f'{self.name}: jet definition: {self.jet_definition.description()}')
        self._log.info(f'{self.name}: jet selector: {self.jet_selector.description()}')

    def process_event(self, input_object, nev=None):
        rv = super(JetAnalysis, self).process_event(input_object, nev=nev)
        self.parts = self.data_source.parts
        if len(self.parts) < 1:
            self._log.debug(f'{self.name}: no particles to analyze len(self.parts)={len(self.parts)}')
            return False
        self.jets = fj.sorted_by_pt(self.jet_selector(self.jet_definition(self.parts)))
        if len(self.jets) < 1:
            return False
        for ij, j in enumerate(self.jets):
            j.set_user_index(ij + self.njet_count)
            self.tn_jet.Fill(self._nev, self.data_source.xsec, self.data_source.ev_weight, self.njet_count + ij, ij, j.pt(), j.eta(), j.phi(), j.m())
        self.njet_count += len(self.jets)
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
            "nev:nj:ij:pt:eta:phi:pt_h2p:xsec:ev_weight:ij_h:pt_h:eta_h:phi_h:pt_ratio:n_matched",
        )
        self.tn_correl_parton_ch = ROOT.TNtuple(
            f"tn_correl_parton_ch_{self.name}",
            "tn_correl_parton_ch",
            "nev:nj:ij:pt:eta:phi:pt_h2p:xsec:ev_weight:ij_h:pt_h:eta_h:phi_h:pt_ratio:n_matched",
        )
        # note the pt cut on jets will be applied to charged jets - the it is somewhat effectively lower for full jets
        self.to_exec = []        
        self.data_source_full = DataSource(self.cfg, name=f"{self.name}_full", source_type=self.data_source.source_type)
        self.data_source_full.select_final(select=True)
        self.data_source_full.select_charged(select=False)
        self.to_exec.append(self.data_source_full)

        self.data_source_charged = DataSource(self.cfg, name=f"{self.name}_charged", source_type=self.data_source.source_type)
        self.data_source_charged.select_charged(select=True)
        self.data_source_charged.select_final(select=True)
        self.to_exec.append(self.data_source_charged)
        self.j_ch_ana = JetAnalysis(self.cfg, name='jet_charged', data_source=self.data_source_charged)
        self.to_exec.append(self.j_ch_ana)

        self.j_full_ana = JetAnalysis(self.cfg, name='jet_full', data_source=self.data_source_full)
        self.to_exec.append(self.j_full_ana)

        if self.cfg.parton_hadron:
            self.data_source_parton = DataSource(self.cfg, name=f"{self.name}_parton", source_type=self.data_source.source_type)
            self.data_source_parton.select_parton(select=True)
            self.to_exec.append(self.data_source_parton)
            self.j_parton_ana = JetAnalysis(self.cfg, name="jet_parton", data_source=self.data_source_parton)
            self.to_exec.append(self.j_parton_ana)

    def process_event(self, input_object, nev=None):
        rv = super(JetChargedFullAnalysis, self).process_event(input_object, nev=nev)
        for e in self.to_exec:
            rv = rv and e.process_event(input_object, nev)
        if not rv:
            return False
        self.tn_events.Fill(self._nev, self.data_source_full.xsec, self.data_source_full.ev_weight,
                            len(self.j_full_ana.jets), len(self.j_ch_ana.jets),
                            len(self.j_full_ana.parts), len(self.j_ch_ana.parts))
        n_matched = {}
        if self.data_source.source_type != 'pythia':
            return rv
        for ij, j in enumerate(self.j_full_ana.jets):
            # this is pythia specific
            _pythia = self.data_source_full.source_object
            sum_charged_in_full = sum([int(_pythia.event[p.user_index()].isCharged()) * p.perp() for p in j.constituents()])
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
                                    self.data_source_full.xsec, self.data_source_full.ev_weight,
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
                    self.tn_correl_parton_full.Fill(self._nev, j.user_index(), ij, j.pt(), j.eta(), j.phi(), sum_charged_in_full,
                                            self.data_source_full.xsec, self.data_source_full.ev_weight,
                                            ijch, jch.pt(), jch.eta(), jch.phi(),
                                            pt_ratio, n_matched[ij])
            n_matched = {}
            for ij, j in enumerate(self.j_parton_ana.jets):
                for ijch, jch in enumerate(self.j_ch_ana.jets):
                    if j.delta_R(jch) > self.cfg.jet_R / 2.0:
                        continue
                    n_matched[ij] = n_matched.get(ij, 0) + 1
                    pt_ratio = jch.pt() / j.pt()
                    self.tn_correl_parton_ch.Fill(self._nev, j.user_index(), ij, j.pt(), j.eta(), j.phi(), sum_charged_in_full,
                                            self.data_source.xsec, self.data_source.ev_weight,
                                            ijch, jch.pt(), jch.eta(), jch.phi(),
                                            pt_ratio, n_matched[ij])

        return rv

class EECAnalysis(AnalysisBase):

    def init(self):
        super(EECAnalysis, self).init()
        self.pt_cuts = [0, 0.15, 1, 2]
        self.tn_eec = {}
        for ptcut in self.pt_cuts:
            self.tn_eec[ptcut] = ROOT.TNtuple(f'tn_eec_{self.name}_ptcut{ptcut}', 'tn_eec', 'nev:xsec:ev_weight:ij:dr:pt1:pt2:eec:ptjet:ptcut')
        if self.jet_analysis is None:
            raise ValueError(f'{self.name}: jet_analysis not set')

    def process_event(self, input_object, nev=None):
        rv = super(EECAnalysis, self).process_event(input_object, nev=nev)
        if not rv:
            return False
        rvs = [self.analyze_eec(ptcut) for ptcut in self.pt_cuts]
        return any(rvs)

    def analyze_eec(self, ptcut):
        # Generate all pairs from parts, excluding pairs of the same element
        # self.pairs = list(itertools.combinations(parts, 2))
        # Generate all pairs from parts, including pairs of the same element
        # _pairs = list(itertools.product(_parts_cut, repeat=2))
        for ij, j in enumerate(self.jet_analysis.jets):
            _parts_cut = [p for p in j.constituents() if p.perp() >= ptcut]
            _pairs = list(itertools.product(_parts_cut, repeat=2))
            log.debug(f'{self.name}: number of pairs: {len(_pairs)} with ptcut: {ptcut}')
            if len(_pairs) < 1:
                continue
            for first, second in _pairs:
                dr = first.delta_R(second)
                eec = first.perp() * second.perp() / pow(j.perp(), 2.)
                self.tn_eec[ptcut].Fill(self._nev, self.jet_analysis.data_source.xsec, self.jet_analysis.data_source.ev_weight, ij, dr, first.perp(), second.perp(), eec, j.perp(), ptcut)
        return True

### PYTHIA CODE FOR D0

def get_D0s_pythia(pythia, skipDstar = True):
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

### HEPMC CODE FOR D0

def get_D0s_skipDstar_hepmc(hepmc_event, noDstar=False):
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
    

def get_D0s_hepmc(hepmc_event, noDstar=False):
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


def get_final_except_hepmc(hepmc_event, parts, charged_only=False):
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

def get_final_hepmc(hepmc_event, charged_only=False):
    fjparts = []
    for i,p in enumerate(hepmc_event.particles):
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


def print_debug_hepmc(Darray, Ddaughters):
    for _p in Darray:
        log.debug(f'D0 is: {_p}')
        for _pp in _p.end_vertex.particles_out:
            log.debug(f'          {_pp}')
        for _pp in Ddaughters:
            log.debug(f'  - check:{_pp}')

