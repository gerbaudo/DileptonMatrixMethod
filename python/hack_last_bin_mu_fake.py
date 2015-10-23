#!/bin/env python
#
# The last high-eta-high-pt bins of the muon fake rate histograms have
# large uncertainties. The down variation can go to zero and make the
# matrix diverge. Patch it by replacing these bin contents with the
# ones from the lower bins (i.e. extrapolate from where we have some
# decent stats).
#

# davide.gerbaudo@gmail.com
# 2015-04-03

import shutil
import ROOT as R

# original_filename ='/tmp/FakeMatrix_Mar_18.root',
# modified_filename ='/tmp/FakeMatrix_Mar_18_hack_mu.root'
# original_filename = '/tmp/gerbaudo/FakeMatrix_Mar_18_syst.root'
# modified_filename = '/tmp/gerbaudo/FakeMatrix_Mar_18_syst_hack_mu.root'
# original_filename = 'data/FakeMatrix_May_10.root'
# modified_filename = 'data/FakeMatrix_May_10_hack_mu.root'
original_filename = 'data/FakeMatrix_Jun_17.root'
modified_filename = 'data/FakeMatrix_Jun_17_hack_mu.root'

shutil.copyfile(original_filename, modified_filename)
input_file = R.TFile(modified_filename, 'update')
h = input_file.Get('mu_fake_rate2d_emuInc')
nx = h.GetNbinsX()
ny = h.GetNbinsY()
for ix in range(1,nx+1):
    for iy in range(1,ny+1):
        print "({},{}) : {:.3f}+/-{:.3f}".format(ix, iy, h.GetBinContent(ix, iy), h.GetBinError(ix, iy))
        if ix==3:
            h.SetBinContent(ix, iy, h.GetBinContent(3,1))
            h.SetBinError  (ix, iy, h.GetBinError(3,1))
print "after hack"
for ix in range(1,nx+1):
    for iy in range(1,ny+1):
        print "({},{}) : {:.3f}+/-{:.3f}".format(ix, iy, h.GetBinContent(ix, iy), h.GetBinError(ix, iy))
h.Write()
