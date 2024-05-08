import itertools
from heppyy.util.logger import Logger
log = Logger()

from heppyy.util.mputils import logbins, filename_safe
import ROOT
import os
import array

class EEChistogram:
	counts_bin = 0
	mult_bin = 1
	npairs = 2
	xsec_bin = 3
	weight_bin = 4
	ev_weight_bin = 5

	def __init__(self, name='eec', binning=18, load=False):
		self.ncorrel = 2
		self.nbins = 18
		self.lbins = logbins(1.e-2, 1., self.nbins)
		if isinstance(binning, int):
			self.nbins = int(binning)
			self.lbins = logbins(1.e-2, 1., self.nbins)
		if isinstance(binning, str):
			self.adopt_binning(binning)
		self.hlist = []
		self.tnlist = []
		self.name = name
		hname = name
		log.debug(f'[d] {name} load is {type(load)}')
		if isinstance(load, ROOT.TFile):
			log.info(f'[i] loading from {load.GetName()}.')
			self.h = load.Get(hname)
			log.info(f'[i] loading {hname+"_pairs"} from {load.GetName()}.')
			self.h_pairs = load.Get(hname+'_pairs')
			log.info(f'[i] loading {hname+"_norm"} from {load.GetName()}.')
			self.h_norm = load.Get(hname+'_norm')
			log.info(f'[i] {self.norms()}.')
			self.hlist = [self.h, self.h_pairs, self.h_norm]
		else:
			log.info(f'[i] creating new histograms... {name}')
			self.h = ROOT.TH1F(hname, hname, self.nbins, self.lbins)
			self.h_pairs = ROOT.TH1F(hname+'_pairs', hname+'_pairs', self.nbins, self.lbins)
			self.h_norm = ROOT.TH1F(hname+'_norm', 'norm', 10, 0., 10.)
			self.hlist = [self.h, self.h_pairs, self.h_norm]
			# add ntuple as x-check - legacy from pythia eech.py
			self.tn = ROOT.TNtuple(hname+'_tn', hname+'tneec', 'n:ptpartcut:RL:w:jetpt:sigmaGen:weight')
			self.tn_eec = ROOT.TNtuple(hname+'_tn_eec', hname+'_tn_eec', 'RL:w:eec:xsec:ew')
			self.tnlist = [self.tn, self.tn_eec]


	def adopt_binning(self, binning):
		nbins = 18
		_lbins = logbins(1.e-2, 1., nbins)
		try:
			_fname = binning.split(':')[0]
			_hname = binning.split(':')[1]
			_lbins = []
			with ROOT.TFile(_fname) as _f:
				htemplate = _f.Get(_hname)
				nbins = htemplate.GetNbinsX()
				for b in range(1, nbins+1):
					_lbins.append(htemplate.GetBinLowEdge(b))
				_lbins.append(htemplate.GetBinLowEdge(nbins) + htemplate.GetBinWidth(b))
				lbins = array.array('f', _lbins)
			log.critical(f'[i] adopting binning {nbins} {lbins}')
		except:
			log.critical(f'[e] unable to get the histogram {binning}')
			return -1
		self.nbins = nbins
		self.lbins = lbins

	def Scale(self, weight, option="width"):
		self.h.Scale(weight, option)
		# do not scale by bin width - just a number of pairs in the bin
		self.h_pairs.Scale(weight)

	def normalize(self, xlow=0.01, xhigh=1.):
		norms = self.norms()
		self.h.Scale(1./norms['ev_weight'])
		_integral_b0 = self.h.FindBin(xlow)
		_integral_b1 = self.h.FindBin(xhigh)
		# _integral = self.h.Integral(_integral_b0, _integral_b1, "width")
		_integral = self.h.Integral(_integral_b0, _integral_b1)
		log.info(f'[w] integral within {xlow}-{xhigh} for {self.name}.{self.h.GetName()} is {_integral}')
		_integral_w = self.h.Integral(_integral_b0, _integral_b1, 'width')
		log.info(f'[w] integral/width within {xlow}-{xhigh} for {self.name}.{self.h.GetName()} is {_integral_w}')
		if _integral == 0:
			_integral = 1.
			log.warning(f'[w] integral is zero within {xlow}-{xhigh} for {self.name}.{self.h.GetName()} - setting it to 1.')
		self.h.Scale(1./_integral, "width")
		self.h_pairs.Scale(1./norms['ncounts']/norms['ev_weight'], "width")

	def fill(self, parts, weight=1., xsec=1., ev_weight=1.):
			self.h_norm.Fill(self.counts_bin, 1.)
			self.h_norm.Fill(self.mult_bin, len(parts))
			self.h_norm.Fill(self.xsec_bin, xsec)
			self.h_norm.Fill(self.weight_bin, weight)
			self.h_norm.Fill(self.ev_weight_bin, ev_weight)
   
			# Generate all pairs from parts, excluding pairs of the same element
			# self.pairs = list(itertools.combinations(parts, 2))
			# Generate all pairs from parts, including pairs of the same element
			self.pairs = list(itertools.product(parts, repeat=2))	
			self.h_norm.Fill(self.npairs, len(self.pairs))
			# Calculate the list of pairs
			self.eec2 = [(first.perp() * second.perp() * weight * ev_weight, first.delta_R(second)) for first, second in self.pairs]
			log.debug(f'number of pairs: {len(self.eec2)}')
			# Fill the histogram
			_ = [self.h.Fill(dr, eec) for eec, dr in self.eec2]
			_ = [self.h_pairs.Fill(dr, ev_weight) for eec, dr in self.eec2]
   
   			# self.tn = ROOT.TNtuple('tn', 'tneec', 'n:ptpartcut:RL:w:jetpt:sigmaGen:weight')
			_ = [self.tn.Fill(2, -1, dr, weight, weight, xsec, ev_weight) for eec, dr in self.eec2]
   			# self.tn_eec = ROOT.TNtuple(hname+'_tn_eec', hname+'_tn_eec', 'RL:w:eec:xsec:ew')
			_ = [self.tn_eec.Fill(dr, weight, eec, xsec, ev_weight) for eec, dr in self.eec2]
   
	def norms(self):
		norms = {}
		norms['ncounts'] = self.h_norm.GetBinContent(self.h_norm.FindBin(self.counts_bin))
		if norms['ncounts'] == 0:
			norms['ncounts'] = 1
			log.critical(f'[w] ncounts is zero for {self.name}.{self.h.GetName()} - setting it to 1.')
		norms['mult'] = self.h_norm.GetBinContent(self.h_norm.FindBin(self.mult_bin))
		norms['npairs'] = self.h_norm.GetBinContent(self.h_norm.FindBin(self.npairs))
		norms['xsec'] = self.h_norm.GetBinContent(self.h_norm.FindBin(self.xsec_bin))
		norms['weight'] = self.h_norm.GetBinContent(self.h_norm.FindBin(self.weight_bin))
		norms['ev_weight'] = self.h_norm.GetBinContent(self.h_norm.FindBin(self.ev_weight_bin))
		norms['mean_xsec'] = norms['xsec'] / norms['ncounts']
		return norms

	def ls(self):
		print(f'  - {self.name}')
		for h in self.hlist:
			_int = h.Integral(h.FindBin(h.GetXaxis().GetXmin()), h.FindBin(h.GetXaxis().GetXmax()), "width")
			_int_uw = h.Integral(h.FindBin(h.GetXaxis().GetXmin()), h.FindBin(h.GetXaxis().GetXmax()))
			print(f'    - {h.GetName()} - entries: {h.GetEntries()} - integral/width: {_int} - integral: {_int_uw}')
		print(f'    - norms: {self.norms()}')

	def Write(self):
		for h in self.hlist:
			h.Write()
		for tn in self.tnlist:
			tn.Write()
  

class EEC2file:
	def __init__(self, output_fname='eec2file_out.root', name='eec2', load=True, use_h_binning=None, args=None):
		self.output_fname = output_fname.replace('.root', '_' + filename_safe(name) + '.root')
		binning = 18
		if use_h_binning:
			binning=use_h_binning
		if load is True:
			self.output_fname = output_fname
			if os.path.exists(self.output_fname) is False:
				log.error(f'[w] file {self.output_fname} does not exist - will create one.')
				# raise FileNotFoundError(f'file {self.output_fname} does not exist.')
				load = False
		if load:
			log.info(f'[i] will load from {self.output_fname}.')
			self.fout = ROOT.TFile(self.output_fname, 'read')
			log.debug(f'[i] will read from {self.fout.GetName()}')
			self.written = True
		else:
			self.fout = ROOT.TFile(self.output_fname, 'recreate')
			self.written = False
		log.debug(f'[i] will write to {self.fout.GetName()}')
		self.fout.cd()
		self.h = {}
		if load:
			self.h['counts'] = self.fout.Get('counts')
			self.h['counts_ew'] = self.fout.Get('counts_w')
			self.h['z_parts_norm'] = self.fout.Get('z_parts_norm')
			self.h['pt_parts_norm'] = self.fout.Get('pt_parts_norm')
			self.h['ev_weights_prof'] = self.fout.Get('ev_weights_prof')
			self.h['xsecs_prof'] = self.fout.Get('xsecs_prof')
		else:
			self.h['counts'] = ROOT.TH1F('counts', 'counts', 200, 0., 200.)
			self.h['counts_ew'] = ROOT.TH1F('counts_w', 'counts_w', 200, 0., 200.)	
			_nbins = 18
			_lbins = logbins(1.e-2, 1., _nbins)
			self.h['z_parts_norm'] = ROOT.TH1F('z_parts_norm', 'z_parts_norm', _nbins, _lbins)
			_lbins = logbins(1.e-2, 100., _nbins)
			self.h['pt_parts_norm'] = ROOT.TH1F('pt_parts_norm', 'pt_parts_norm', _nbins, _lbins)

		_fload = load
		if load is True and self.fout:
			_fload = self.fout
		self.h['eec2_pt_0.0'] = EEChistogram(name='eec2_pt0', binning=binning, load=_fload)
		self.h['eec2_pt_0.15'] = EEChistogram(name='eec2_pt0.15', binning=binning, load=_fload)
		self.h['eec2_pt_1.0'] = EEChistogram(name='eec2_pt1', binning=binning, load=_fload)

		self.h['eec2_pt_0.0_ew'] = EEChistogram(name='eec2_pt0_ew', binning=binning, load=_fload)
		self.h['eec2_pt_0.15_ew'] = EEChistogram(name='eec2_pt0.15_ew', binning=binning, load=_fload)
		self.h['eec2_pt_1.0_ew'] = EEChistogram(name='eec2_pt1_ew', binning=binning, load=_fload)
  
		self.ev_weights = []
		self.xsecs = []
  
	def fill(self, parts, pTweight, xsection, event_weight=1.):
		self.ev_weights.append(event_weight)
		self.xsecs.append(xsection)
  
		self.h['counts'].Fill(pTweight)
		self.h['counts_ew'].Fill(pTweight, event_weight)
 
		# fill the z for the object (jet if jet passed)
		_ = [self.h['z_parts_norm'].Fill(p.perp() / pTweight) for p in parts]
		_ = [self.h['pt_parts_norm'].Fill(p.perp()) for p in parts]
 
		_pTweight2 = pTweight * pTweight
		_parts = parts
		self.h['eec2_pt_0.0'].fill(_parts, 		weight=1./_pTweight2, xsec=xsection, ev_weight=1.)
		self.h['eec2_pt_0.0_ew'].fill(_parts, 	weight=1./_pTweight2, xsec=xsection, ev_weight=event_weight)

		_parts = [p for p in parts if p.perp() > 0.15]
		self.h['eec2_pt_0.15'].fill(_parts, 	weight=1./_pTweight2, xsec=xsection, ev_weight=1.)
		self.h['eec2_pt_0.15_ew'].fill(_parts, 	weight=1./_pTweight2, xsec=xsection, ev_weight=event_weight)

		_parts = [p for p in parts if p.perp() > 1.0]
		self.h['eec2_pt_1.0'].fill(_parts, 		weight=1./_pTweight2, xsec=xsection, ev_weight=1.)
		self.h['eec2_pt_1.0_ew'].fill(_parts, 	weight=1./_pTweight2, xsec=xsection, ev_weight=event_weight)

	def normalize(self, xlow=0.01, xhigh=1.):
		norm_counts = self.h['counts'].Integral()
		norms = self.h['eec2_pt_0.0'].norms()
		log.debug(f'norms: {norms["ncounts"]} - {norm_counts}')  
		for h in self.h.values():
			if isinstance(h, EEChistogram):
				h.normalize(xlow=xlow, xhigh=xhigh)
			# normalize the histograms
			if isinstance(h, ROOT.TH1F):
				if 'norm' in h.GetName():
					if norm_counts > 0:
						if 'bw' in h.GetName():
							h.Scale(1./norm_counts, "width")
						else:
							h.Sumw2()
							h.Scale(1./norm_counts)

	def norms(self):
		norms = self.h['eec2_pt_0.0'].norms()
		return norms

	def ls(self):
		print(f'- {self.output_fname}')
		for h in self.h.values():
			if isinstance(h, EEChistogram):
				h.ls()
			else:
				print(f'  - {h.GetName()} - {h.GetEntries()} entries - {h.Integral()} integral.')

	def Write(self):
		if self.fout:
			self.fout.cd()
			for h in self.h.values():
				h.Write()
			# write the weights and xsecs
			p_weights = ROOT.TProfile('ev_weights_prof', 'ev_weights_prof', len(self.ev_weights), 0., len(self.ev_weights))
			_ = [p_weights.Fill(i, w) for i, w in enumerate(self.ev_weights)]
			p_xsecs = ROOT.TProfile('xsecs_prof', 'xsecs_prof', len(self.xsecs), 0., len(self.xsecs))
			_ = [p_xsecs.Fill(i, x) for i, x in enumerate(self.xsecs)]
			self.fout.Write()
			self.fout.Close()
			log.critical(f'[i] written {self.fout.GetName()}.')
		else:
			log.debug(f'[i] nothing to write. self.fout is {self.fout}.')
		self.written = True

	def SaveAs(self, fname):
		_output_fname = fname
		_fout = ROOT.TFile(_output_fname, 'recreate')
		_fout.cd()
		for h in self.h.values():
			if h is not None:				
				h.Write()
		_fout.Write()
		log.critical(f'[i] written {_fout.GetName()}.')
		_fout.Close()
  
	def Close(self):
		if self.fout:
			self.fout.Close()

	def __del__(self):
		if self.written is False:
			log.critical(f'[w] will try to write output file...')
			self.Write()
		if self.fout:
			self.fout.Close()
  
if __name__ == '__main__':
	log.set_level('DEBUG')
	halt = EEC2file(output_fname='example_eec.root', name='eec2', args=None)
	# then for each jet in the event, fill the histogram
	# for j in jets:
	# 	_pTweight = j.perp()
	# 	_parts = j.constituents()
	# 	_parts = [p for p in _parts if p.perp() > 0.15]
	#   halt.fill(_parts, _pTweight, _event_weight)
	del halt