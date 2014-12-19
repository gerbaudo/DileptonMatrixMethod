#!/bin/env python

# Perform a closure test for a given input fake matrix
#
# Given a fake matrix, retrieve the fake and real efficiencies for
# electron and muon. Then generate dilepton pairs that are TT/TL/LT/LL
# according to these efficiencies. Feed these dilepton pairs to the
# DileptonMatrixMethod class, and verify that you get the number of
# fake that were generated.
#
# davide.gerbaudo@gmail.com
# 2013-12-18

# for now always generate emu
# for now only consider the 1D parametrization vs pt

import collections
import pprint
import sys
from random import random, randint

import ROOT as r
r.gROOT.SetBatch(1)
from print_available_regions import get_histonames, classify_by_region
from compare_matrices import getBinIndices, binContents, formatted_histo_contents, bc_format
from run_on_txt_input import load_packages

fraction_of_real = 0.25 # this depends on how much of a fake-enriched region you're simulating

class Efficiency:
    def __init__(self, n_pt_bins):
        self.num = n_pt_bins*[0.0]
        self.den = n_pt_bins*[0.0]
    def count(self, pt_bin, is_tight):
        self.den[pt_bin] += 1.0
        self.num[pt_bin] += (1.0 if is_tight else 0.0)
    def values(self):
        return [(n/d if d else 0.0) for n,d in zip(self.num, self.den)]
    def values_as_str(self):
        return ' '.join(bc_format(v, 4) for v in self.values())

def main():
    if len(sys.argv)<2:
        print "Usage:\n%s matrix.root"%sys.argv[0]
        sys.exit(1)
    matrix_name = sys.argv[1]
    input_file = r.TFile.Open(matrix_name)
    if not input_file:
        print "invalid matrix"
        return
    print "reading efficiencies from {0}".format(input_file.GetName())
    histonames_per_region = classify_by_region(get_histonames(input_file))['complete']
    print "available regions : {0}".format(' '.join(histonames_per_region.keys()))

    regions = histonames_per_region.keys()
    regions2d = [region for region in regions if '2d' in region]
    regions1d = [region for region in regions if region not in regions2d]
    histos_per_region = dict((region, [input_file.Get(hn) for hn in hns])
                             for region, hns in histonames_per_region.iteritems())
    region = 'emuInc' # for now choose this one
    if region not in regions1d:
        print "region {0} is not available".format(region)
        return
    def find_by_name(objects=[], name=''):
        obj = None
        for o in objects:
            if name==o.GetName():
                obj = o
                break
        return obj

    h_eff_el_real = find_by_name(histos_per_region[region], 'el_real_eff_' +region)
    h_eff_el_fake = find_by_name(histos_per_region[region], 'el_fake_rate_'+region)
    h_eff_mu_real = find_by_name(histos_per_region[region], 'mu_real_eff_' +region)
    h_eff_mu_fake = find_by_name(histos_per_region[region], 'mu_fake_rate_'+region)
    assert(1==len(set(h.GetNbinsX() for h in [h_eff_el_real, h_eff_el_fake, h_eff_mu_real, h_eff_mu_fake]))),"more than one binning"
    pt_bin_indices = getBinIndices(h_eff_el_real)
    n_pt_bins = len(pt_bin_indices)
    pt_min = h_eff_el_real.GetBinLowEdge(pt_bin_indices[0])
    pt_max = h_eff_el_real.GetBinLowEdge(pt_bin_indices[-1])+h_eff_el_real.GetBinWidth(pt_bin_indices[-1])
    pt_bin_edges = [h_eff_el_real.GetBinLowEdge(b) for b in pt_bin_indices] + [pt_max]
    
    eff_el_real = binContents(h_eff_el_real)
    eff_el_fake = binContents(h_eff_el_fake)
    eff_mu_real = binContents(h_eff_mu_real)
    eff_mu_fake = binContents(h_eff_mu_fake)


    load_packages()
    matrix = r.susy.fake.DileptonMatrixMethod()
    regions = r.std.vector('std::string')()
    regions.push_back(region)
    matrix.configure1d(matrix_name, regions)
    systematic = r.susy.fake.Systematic.SYS_NOM
    region_index = matrix.getIndexRegion(region)
    l0 = r.susy.fake.Lepton(False, False, 10.0, 1.0)
    l1 = r.susy.fake.Lepton(False, False, 10.0, 1.0)

    


    num_events_to_generate = int(1e6)
    num_events_to_print = 10
    num_TT = 0
    num_fake_TT = 0
    num_FF_TT, num_FR_TT, num_RF_TT, num_RR_TT = 0, 0, 0, 0
    est_fake_TT = 0.0
    eff_el_real_meas = Efficiency(n_pt_bins)
    eff_el_fake_meas = Efficiency(n_pt_bins)
    eff_mu_real_meas = Efficiency(n_pt_bins)
    eff_mu_fake_meas = Efficiency(n_pt_bins)

    def random_npt():
        "normalized pt, between 0 and 1"
        return pt_values[randint(0, n_pt_bins)]
    def random_is_real():
        return random()<fraction_of_real
    for i in xrange(num_events_to_generate):
        el_npt = random()
        mu_npt = random()
        el_pt = pt_min + el_npt*pt_max
        mu_pt = pt_min + mu_npt*pt_max
        el_ipt = int(el_npt*n_pt_bins)
        mu_ipt = int(mu_npt*n_pt_bins)
        el_is_real = random_is_real()
        mu_is_real = random_is_real()
        el_is_fake = not el_is_real
        mu_is_fake = not mu_is_real
        el_eff = eff_el_real[el_ipt] if el_is_real else eff_el_fake[el_ipt]
        mu_eff = eff_mu_real[mu_ipt] if mu_is_real else eff_mu_fake[mu_ipt]
        el_is_tight = random()<el_eff
        mu_is_tight = random()<mu_eff
        is_fake_event = (el_is_fake or mu_is_fake)

        if el_is_real : eff_el_real_meas.count(el_ipt, el_is_tight)
        else          : eff_el_fake_meas.count(el_ipt, el_is_tight)
        if mu_is_real : eff_mu_real_meas.count(mu_ipt, mu_is_tight)
        else          : eff_mu_fake_meas.count(mu_ipt, mu_is_tight)

        l0 = l0.isEl(True ).isTight(el_is_tight).pt(el_pt)#.eta(0.1)
        l1 = l1.isEl(False).isTight(mu_is_tight).pt(mu_pt)#.eta(0.1)
        l0, l1 = ((l0, l1) if el_pt>mu_pt else (l1, l0))
        weight = matrix.getTotalFake(l0, l1, region_index, systematic)
        if el_is_tight and mu_is_tight:
            num_TT += 1            
            num_fake_TT += (1 if is_fake_event else 0)
            est_fake_TT += weight
            num_FF_TT += (1 if (el_is_fake and mu_is_fake) else 0)
            num_FR_TT += (1 if (el_is_fake and mu_is_real) else 0)
            num_RF_TT += (1 if (el_is_real and mu_is_fake) else 0)
            num_RR_TT += (1 if (el_is_real and mu_is_real) else 0)

        if i<num_events_to_print:
            print "gen el({0}) mu({1}) -> l0 {2}, l2 {3}, w = {4}".format('R' if el_is_real else 'F',
                                                                          'R' if mu_is_real else 'F',
                                                                          l0.str(), l1.str(), weight)
    print "summary"
    print "generated: {0} events".format(num_events_to_generate)
    print "of which {0} were TT".format(num_TT)
    print "of which {0} were due to fake".format(num_fake_TT)
    print 'num_FF_TT ',num_FF_TT
    print 'num_FR_TT ',num_FR_TT
    print 'num_RF_TT ',num_RF_TT
    print 'num_RR_TT ',num_RR_TT
    print "estimated number of fakes: {0}".format(est_fake_TT)
    print "efficiencies retrieved:"
    for h in [h_eff_el_real, h_eff_el_fake, h_eff_mu_real, h_eff_mu_fake]:
        print "{0} : {1}".format(h.GetName(), formatted_histo_contents(h))
    print "efficiencies measured:"
    for h in ['eff_el_real_meas', 'eff_el_fake_meas', 'eff_mu_real_meas', 'eff_mu_fake_meas']:
        print "{0} : {1}".format(h, eval(h).values_as_str())

#___________________________________________________________

main()
