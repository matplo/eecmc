import itertools
from heppyy.util.logger import Logger
log = Logger()

from heppyy.util.mputils import logbins, filename_safe
import ROOT

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

		self.h['eec2_pt_0.0_ew'] = EEChistogram(name='eec2_pt0_ew')
		self.h['eec2_pt_0.15_ew'] = EEChistogram(name='eec2_pt0.15_ew')
		self.h['eec2_pt_1.0_ew'] = EEChistogram(name='eec2_pt1_ew')
  
		self.weights = []
		self.xsecs = []

	def fill(self, parts, pTweight, xsection, event_weight=1.):
		self.weights.append(event_weight)
		self.xsecs.append(xsection)
  
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

	def __del__(self):
		if self.fout:
			fname = self.fout.GetName()
			self.fout.cd()
			norm = self.h['counts'].Integral()
			norm = len(self.xsecs) * 1.
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
			# write the weights and xsecs
			p_weights = ROOT.TProfile('weights', 'weights', len(self.weights), 0., len(self.weights))
			_ = [p_weights.Fill(i, w) for i, w in enumerate(self.weights)]
			p_xsecs = ROOT.TProfile('xsecs', 'xsecs', len(self.xsecs), 0., len(self.xsecs))
			_ = [p_xsecs.Fill(i, x) for i, x in enumerate(self.xsecs)]
			self.fout.Write()
			self.fout.Close()
			log.debug(f'[i] written {fname}.')
		else:
			log.debug(f'[i] nothing to write. self.fout is {self.fout}.')

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