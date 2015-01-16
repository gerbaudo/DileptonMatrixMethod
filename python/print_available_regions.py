#!/bin/env python

# Script to print out the list of fake histograms for each signal region
#
# davide.gerbaudo@gmail.com
# July 2014

import collections
import pprint
import sys
import ROOT as r
r.gROOT.SetBatch(1)

def main():
    if len(sys.argv)<2:
        print "Usage:\n%s filename_matrix.root"%sys.argv[0]
        sys.exit(1)
    filename = sys.argv[1]
    input_file = r.TFile.Open(filename)
    if input_file:
        histonames = get_histonames(input_file)
        histos_by_region = classify_by_region(histonames)        
        pprint.pprint(histos_by_region)
#___________________________________________________________
                
def get_histonames(input_file):
    keys = input_file.GetListOfKeys()
    return list(set([k.GetName() for k in keys if isTH(k.GetClassName())]))

def isTH(classname):
    if not hasattr(isTH, 'th1') : isTH.th1 = r.TH1.Class() # cache function attr
    return r.TClass(classname).InheritsFrom(isTH.th1)
def isTH1(classname):
    if not hasattr(isTH1, 'th1') : isTH1.th1 = r.TH1.Class()
    if not hasattr(isTH1, 'th2') : isTH1.th2 = r.TH2.Class()
    cl = r.TClass(classname)
    return cl.InheritsFrom(isTH1.th1) and not cl.InheritsFrom(isTH1.th2)
def isTH2(classname):
    if not hasattr(isTH2, 'th2') : isTH2.th2 = r.TH2.Class()
    cl = r.TClass(classname)
    return cl.InheritsFrom(isTH2.th2)

def classify_by_region(histonames=[]):
    def extract_region(hn):
        "see DileptonMatrixMethod::generateNominalHistoname"
        return (hn
                .replace('el_real_eff_','')
                .replace('el_fake_rate_', '')
                .replace('mu_real_eff_', '')
                .replace('mu_fake_rate_', '')
                .replace('el_real_eff2d_', '')
                .replace('el_fake_rate2d_', '')
                .replace('mu_real_eff2d_', '')
                .replace('mu_fake_rate2d_', ''))
    histos_by_region = collections.defaultdict(list)
    for hn in histonames:
        histos_by_region[extract_region(hn)].append(hn)
    histos_by_cat = {'complete' : dict(),
                     'incomplete' : dict()}
    def region_is_complete(histos_this_region=[]):
        # expect either 4 (1D or 2D) or 8 histos (both)
        return len(histos_this_region) in [4, 8]
    for r, histos in histos_by_region.iteritems():
        if region_is_complete(histos):
            histos_by_cat['complete'][r] = histos
        else:
            histos_by_cat['incomplete'][r] = histos
    return histos_by_cat

if __name__=='__main__':
    main()
