#!/bin/env python

# Script to print out the rates from two matrix files
#
# Only the regions that are available in both files are considered.
# If only one file is provided (twice), a summary table with all
# values is printed out.
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
    if (filename1 == filename2) and input_file1:
        print_summary_table(input_file1)
    elif input_file1 and input_file2:
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

def bc_format(bc, precision):
    # return ('{:<'+str(precision+3)+'.'+str(precision)+'}').format(bc)
    return ('{:.'+str(precision)+'}').format(bc).ljust(precision+len('0.0'))
def ratio_format(a, b, precision):
    width = precision+len('0.0')
    return (('{:.'+str(precision)+'}').format(a/b) if b else '--').ljust(width)
def formatted_histo_contents(h, precision=4):
    return ' '.join(bc_format(bc, precision) for bc in binContents(h))
def formatted_histos_ratio(h_num, h_den, precision=4):
    return ' '.join(ratio_format(n, d, precision) for n,d in zip(binContents(h_num), binContents(h_den)))

def print_histos_comparison(common_histos={'region':[]}, input_file1=None, input_file2=None):
    for region, histos in common_histos.iteritems():
        print
        print '--- '+region+' ---'
        for hn in histos:
            h1 = input_file1.Get(hn)
            h2 = input_file2.Get(hn)
            print
            print '--  '+hn+'  --'
            print "1)   {0}".format(formatted_histo_contents(h1) if h1 else '--')
            print "2)   {0}".format(formatted_histo_contents(h2) if h2 else '--')
            print "2/1) {0}".format(formatted_histos_ratio(h2, h1))

def print_summary_table(input_file):
    "Loop over all histograms in the file, and summarized them in a table"
    print "Summary for {0}".format(input_file.GetName())
    histonames_per_region = classify_by_region(get_histonames(input_file))['complete']
    regions = histonames_per_region.keys()
    regions2d = [r for r in regions if '2d' in r]
    regions1d = [r for r in regions if r not in regions2d]
    histos_per_region = dict((r, [input_file.Get(hn) for hn in hns])
                             for r, hns in histonames_per_region.iteritems())
    # the line lenght is quite different for 1D and 2D histos, so print them separately
    def is_2d(h) : return '2d' in h.GetName()
    def is_1d(h) : return not is_2d(h)
    histos1d = dict((r, filter(is_1d, histos_per_region[r])) for r in regions)
    histos2d = dict((r, filter(is_2d, histos_per_region[r])) for r in regions)
    if histos1d:
        print "1D"
        print_histos_contents_table(histos1d)
    if histos2d:
        print "2D"
        print_histos_contents_table(histos2d)

def print_histos_contents_table(histos_per_region={}):
    "Print a summary with one line per region"
    histonames_per_region = dict((r, [h.GetName() for h in histos])
                                 for r, histos in histos_per_region.iteritems())
    histoname_headers = sorted(list(set(h.replace('_'+r, '')
                                        for r, hs in histonames_per_region.iteritems()
                                        for h in hs)))
    if not histoname_headers:
        print "no histograms"
        return
    def h_matches_this_col(histo, col_header, region):
        return col_header==histo.GetName().replace('_'+region, '')
    # sort the histos in the dict according to the header
    histos_per_region = dict((r, [[h for h in histos_per_region[r]
                                   if h_matches_this_col(h, hh, r)][0]
                                  for hh in histoname_headers])
                             for r in histos_per_region.keys())
    precision = 4
    max_width_first_col = max(len(f) for f in histoname_headers)
    max_width_other_col = max(len(f) for f in
                              [formatted_histo_contents(h, precision)
                               for hns in histos_per_region.values() for h in hns])
    header = '  '.join(['region'.ljust(max_width_first_col)] +
                       [h.ljust(max_width_other_col) for h in histoname_headers])
    lines =  '\n'.join(['  '.join([r.ljust(max_width_first_col)] +
                                  [formatted_histo_contents(h).ljust(max_width_other_col)
                                   for h in histos_per_region[r]])
                        for r in sorted(histos_per_region.keys())])
    print header
    print lines


if __name__=='__main__':
    main()
