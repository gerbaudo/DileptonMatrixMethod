#!/bin/env python

# Script to print out the rates from two matrix files
#
# Only the regions that are available in both files are considered.
#
# davide.gerbaudo@gmail.com
# Nov 2014

import collections
import pprint
import sys
import ROOT as r
r.gROOT.SetBatch(1)
from print_available_regions import get_histonames, classify_by_region

def main():
    if len(sys.argv)<3:
        print "Usage:\n%s matrix1.root matrix2.root"%sys.argv[0]
        sys.exit(1)
    filename1 = sys.argv[1]
    filename2 = sys.argv[2]
    input_file1 = r.TFile.Open(filename1)
    input_file2 = r.TFile.Open(filename2)
    if input_file1 and input_file2:
        print "inputs:"
        print "1) {0}".format(input_file1.GetName())
        print "2) {0}".format(input_file2.GetName())
        histos1 = classify_by_region(get_histonames(input_file1))['complete']
        histos2 = classify_by_region(get_histonames(input_file2))['complete']
        common_regions = list(set(histos1.keys()).intersection(histos2.keys()))
        common_histos = dict((r, sorted(list(set(histos1[r]).intersection(histos2[r]))))
                             for r in common_regions)
        print_histos_comparison(common_histos, input_file1, input_file2)
#___________________________________________________________
def getBinIndices(h) :
    "Return a list of the internal indices used by TH1/TH2/TH3; see TH1::GetBin for info on internal mapping"
    cname = h.Class().GetName()
    if   cname.startswith('TH1') :
        return [h.GetBin(i)
                for i in range(1, 1+h.GetNbinsX())]
    elif cname.startswith('TH2') :
        return [h.GetBin(i, j)
                for i in range(1, 1+h.GetNbinsX())
                for j in range(1, 1+h.GetNbinsY())]
    elif cname.startswith('TH3') :
        return [h.GetBin(i, j, k)
                for i in range(1, 1+h.GetNbinsX())
                for j in range(1, 1+h.GetNbinsY())
                for k in range(1, 1+h.GetNbinsZ())]
    else : return []

def binContents(h) :
    return [h.GetBinContent(b) for b in getBinIndices(h)]

def print_histos_comparison(common_histos={'region':[]}, input_file1=None, input_file2=None):
    def formatted_histo_contents(h, precision=4):
        def bc_format(bc):
            # return ('{:<'+str(precision+3)+'.'+str(precision)+'}').format(bc)
            return ('{:.'+str(precision)+'}').format(bc).ljust(precision+len('0.0'))
        return ' '.join(bc_format(bc) for bc in binContents(h))
    for region, histos in common_histos.iteritems():
        print
        print '--- '+region+' ---'
        for hn in histos:
            h1 = input_file1.Get(hn)
            h2 = input_file2.Get(hn)
            print
            print '--  '+hn+'  --'
            print "1) {0}".format(formatted_histo_contents(h1) if h1 else '--')
            print "2) {0}".format(formatted_histo_contents(h2) if h2 else '--')

if __name__=='__main__':
    main()
