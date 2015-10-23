#!/bin/env python

# Script to plot the emuInc efficiencies
#
# davide.gerbaudo@gmail.com
# July 2014

import collections
import pprint
import sys
import ROOT as r
r.gROOT.SetBatch(1)
r.gROOT.SetStyle('Plain')
r.gStyle.SetPadTickX(1)
r.gStyle.SetPadTickY(1)

from print_available_regions import get_histonames, classify_by_region

def main():
    if len(sys.argv)<2:
        print "Usage:\n%s filename_matrix.root"%sys.argv[0]
        sys.exit(1)
    filename = sys.argv[1]
    input_file = r.TFile.Open(filename)
    can = r.TCanvas('c','')
    if input_file:
        histonames = get_histonames(input_file)
        histos_by_region = classify_by_region(histonames)['complete']
        for region, histonames in histos_by_region.iteritems():
            print '-- ',region,' --'
            histonames = [hn for hn in histonames if '2d' not in hn]
            print histonames
            for hn in histonames:
                h = input_file.Get(hn)
                h.SetTitle(titles[hn])
                h.GetXaxis().SetTitle('p_{T} [GeV]')
                h.GetYaxis().SetTitle('#varepsilon(T|L)')
                h.GetYaxis().SetRangeUser(0.0, 1.1)
                if h.GetMinimum() > 0.0:
                    h.GetYaxis().SetRangeUser(0.0, 1.1)
                h.SetStats(False)
                can.cd()
                can.Clear()
                h.Draw()
                can.Update()
                can.SaveAs(h.GetName()+'.eps')
                can.SaveAs(h.GetName()+'.png')
            # now plot them together
            can.cd()
            can.Clear()
            h_real = input_file.Get('el_real_eff_emuInc')
            h_fake = input_file.Get('el_fake_rate_emuInc')
            h_real.SetTitle(h_real.GetTitle().replace(' real', ''))
            h_real.Draw()
            h_fake.SetMarkerColor(r.kRed)
            h_fake.Draw('same')
            can.SaveAs('el_emuInc.png')

            can.cd()
            can.Clear()
            h_real = input_file.Get('mu_real_eff_emuInc')
            h_fake = input_file.Get('mu_fake_rate_emuInc')
            h_real.SetTitle(h_real.GetTitle().replace(' real', ''))
            h_real.Draw()
            h_fake.SetMarkerColor(r.kRed)
            h_fake.Draw('same')
            can.SaveAs('mu_emuInc.png')
            
            
titles = {
    'el_fake_rate_emuInc' : ' #varepsilon(T|L) fake el',
    'mu_real_eff_emuInc' : ' #varepsilon(T|L) real mu',
    'mu_fake_rate_emuInc' : ' #varepsilon(T|L) fake mu',
    'el_real_eff_emuInc' : ' #varepsilon(T|L) real el',
    }
#___________________________________________________________
if __name__=='__main__':
    main()
