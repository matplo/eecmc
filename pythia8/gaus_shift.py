#!/usr/bin/env python

import ROOT

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, 'analysis'))
from analysis import SingleRootFile

rf = SingleRootFile(fname='gaus_shift.root')
rf.root_file.cd()

# code demonstrating ratio of two gaussian functions where one is shifted by epsilon - parameters otherwise the same

# create two gaussian functions
g1 = ROOT.TF1("g1","TMath::Gaus(x,0,1)",-5,5)
g2 = ROOT.TF1("g2","TMath::Gaus(x,0.1,1)",-5,5)

# create a canvas
c = ROOT.TCanvas("c","c",800,600)
c.Divide(1,2)

c.cd(1)
# draw the first gaussian
g1.Draw()
g1.SetLineColor(ROOT.kRed)
g1.SetLineStyle(2)

# draw the second gaussian
g2.Draw("same")
g2.SetLineColor(ROOT.kBlue)
g2.SetLineStyle(2)

# create a new function that is the ratio of the two gaussians
ratio = ROOT.TF1("ratio","g2/g1",-5,5)
ratio.SetLineColor(ROOT.kGreen)
ratio.SetLineStyle(2)

c.cd(2)
# draw the ratio
ratio.Draw("")

# fit the ratio with pol(2) function
hratio = ROOT.TH1F("hratio","hratio",100,-5,5)
ntrials = 1000000
hratio.FillRandom("ratio", ntrials)
hratio.Scale(1./ntrials)
hratio.Draw("same")
fpol2 = ROOT.TF1("fpol2","pol2",-5,5)
hratio.Fit("fpol2")

for i in range(fpol2.GetNpar()):
	print(i, fpol2.GetParameter(i))

# draw the canvas
c.Update()
c.SaveAs("ratio.png")

for o in [g1,g2,ratio,hratio,fpol2]:
  rf.root_file.cd()
  o.Write()
rf.close()

# wait for input to keep the GUI (which lives on a ROOT event dispatcher) alive
if __name__ == "__main__":
		rep = ''
		while not rep in [ 'q', 'Q' ]:
				rep = input( 'enter "q" to quit: ' )
				if 1 < len(rep):
						rep = rep[0]

