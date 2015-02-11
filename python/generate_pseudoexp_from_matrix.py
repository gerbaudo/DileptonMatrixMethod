#!/bin/env python

# Perform a closure test for a given input fake matrix
#
# Given a fake matrix, retrieve the fake and real efficiencies for
# electron and muon. Then generate dilepton pairs that are TT/TL/LT/LL
# according to these efficiencies. Feed these dilepton pairs to the
# DileptonMatrixMethod class, and verify that you get the number of
# fake that were generated (modulo statistical fluctuations and modulo
# the limitations of the matrix method, discussed in hep-ph/1407.5624).
#
# davide.gerbaudo@gmail.com
# 2013-12-18

# For now only consider the 1D parametrization vs pt; the pt-eta
# parametrization is the same, provided that the efficiency lookup is
# done correctly.

import collections
import pprint
import sys
from random import random, randint

import ROOT as r
r.gROOT.SetBatch(1)
from print_available_regions import get_histonames, classify_by_region
from compare_matrices import getBinIndices, binContents, formatted_histo_contents, bc_format
from run_on_txt_input import load_packages

# this depends on how much of a fake-enriched region you're simulating
fraction_of_real_l0 = 0.95
fraction_of_real_l1 = 0.25

class Efficiency:
    "An object to keep track of the tight/loose efficiency"
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
    if len(sys.argv)<3:
        print "Usage:\n%s matrix.root 10000"%sys.argv[0]
        sys.exit(1)
    matrix_name = sys.argv[1]
    num_TT_events_to_generate = int(sys.argv[2])
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
    def random_is_real(fraction_of_real):
        return random()<fraction_of_real
    def random_is_tight(loose_to_tight_efficiency):
        return random()<loose_to_tight_efficiency


    print "\n"
    print "Test of single-category matrix:"
    # first test: one lepton, one real eff and one fake eff (i.e. one
    # pt bin) use electron efficiencies, but any combination of el/mu
    # would be fine to check the getTotalFake function.
    r0 = r1 = eff_el_real[0]
    f0 = f1 = eff_el_fake[0]
    num_RR, num_RF, num_FR, num_FF = 0, 0, 0, 0
    num_TT, num_TL, num_LT, num_LL = 0, 0, 0, 0
    num_FF_LL, num_FR_LL, num_RF_LL, num_RR_LL = 0, 0, 0, 0
    num_FF_TT, num_FR_TT, num_RF_TT, num_RR_TT = 0, 0, 0, 0
    # measured efficiencies
    m_r0, m_r1 = Efficiency(1), Efficiency(1)
    m_f0, m_f1 = Efficiency(1), Efficiency(1)
    num_events_generated = 0
    while num_TT < num_TT_events_to_generate:
        num_events_generated += 1
        l0_is_real = random_is_real(fraction_of_real_l0)
        l1_is_real = random_is_real(fraction_of_real_l1)
        l0_is_fake = not l0_is_real
        l1_is_fake = not l1_is_real
        l0_is_tight = random_is_tight(r0 if l0_is_real else f0)
        l1_is_tight = random_is_tight(r1 if l1_is_real else f1)
        l0_is_loose = not l0_is_tight
        l1_is_loose = not l1_is_tight
        RR = (l0_is_real and l1_is_real)
        RF = (l0_is_real and l1_is_fake)
        FR = (l0_is_fake and l1_is_real)
        FF = (l0_is_fake and l1_is_fake)
        TT = (l0_is_tight and l1_is_tight)
        TL = (l0_is_tight and l1_is_loose)
        LT = (l0_is_loose and l1_is_tight)
        LL = (l0_is_loose and l1_is_loose)
        num_RR += 1 if RR else 0
        num_RF += 1 if RF else 0
        num_FR += 1 if FR else 0
        num_FF += 1 if FF else 0
        num_TT += 1 if TT else 0
        num_TL += 1 if TL else 0
        num_LT += 1 if LT else 0
        num_LL += 1 if LL else 0
        num_FF_TT += 1 if (FF and TT) else 0
        num_FR_TT += 1 if (FR and TT) else 0
        num_RF_TT += 1 if (RF and TT) else 0
        num_RR_TT += 1 if (RR and TT) else 0
        num_FF_LL += 1 if FF else 0
        num_FR_LL += 1 if FR else 0
        num_RF_LL += 1 if RF else 0
        num_RR_LL += 1 if RR else 0
        if l0_is_real : m_r0.count(0, l0_is_tight)
        else          : m_f0.count(0, l0_is_tight)
        if l1_is_real : m_r1.count(0, l1_is_tight)
        else          : m_f1.count(0, l1_is_tight)

    print "Summary of the test:"
    print "generated: {0} events".format(num_events_generated)
    print "of which {0} were TT".format(num_TT)
    print "of which {0} were due to fake ({1} FF + {2} FR + {3} RF)".format(num_FF_TT + num_FR_TT + num_RF_TT,
                                                                            num_FF_TT, num_FR_TT, num_RF_TT)
    print "input efficiencies    : r0 {0:.0%}, r1 {1:.0%}, f0 {2:.0%}, f1 {3:.0%}".format(r0, r1, f0, f1)
    print "measured efficiencies : r0 {0:.0%}, r1 {1:.0%}, f0 {2:.0%}, f1 {3:.0%}".format(m_r0.values()[0], m_r1.values()[0],
                                                                                          m_f0.values()[0], m_f1.values()[0])
    m_eff = matrix_from_efficiencies(r0, r1, f0, f1)
    m_LL = matrix_LL_events_compositions(num_RR_LL, num_RF_LL, num_FR_LL, num_FF_LL)
    # m_res = r.TMatrixD(m_input, r.TMatrixD.kMult, m_LL)
    m_selected = m_eff * m_LL
    print "n_TT from matrix multiplication (using input    efficiencies) : {0:.3f}".format(m_selected[0][0])
    m_meas_eff = matrix_from_efficiencies(m_r0.values()[0], m_r1.values()[0], m_f0.values()[0], m_f1.values()[0])
    m_selected = m_meas_eff * m_LL
    print "n_TT from matrix multiplication (using measured efficiencies) : {0:.3f}".format(m_selected[0][0])

    est_fake_TT = matrix.getTotalFake(num_TT, num_TL, num_LT, num_LL, r0, r1, f0, f1)
    m_selected = matrix_selected_events(num_TT, num_TL, num_LT, num_LL)
    est_fake_TT = n_fake_from_efficiencies_and_selected(r0, r1, f0, f1, m_selected)
    print "true number of fakes:  {0}".format(num_FF_TT + num_FR_TT + num_RF_TT)
    print "from matrix tool:      {0:.3f}".format(est_fake_TT)
    print "from matrix inversion: {0:.3f}".format(est_fake_TT)

    print "\n"
    print "Test of multiple-category matrix:"
    print "Now doing the same thing, but generate el+mu pairs, and span several pt bins"
    tot_fake_weight = 0.0
    num_RR, num_RF, num_FR, num_FF = 0, 0, 0, 0
    num_TT, num_TL, num_LT, num_LL = 0, 0, 0, 0
    num_FF_LL, num_FR_LL, num_RF_LL, num_RR_LL = 0, 0, 0, 0
    num_FF_TT, num_FR_TT, num_RF_TT, num_RR_TT = 0, 0, 0, 0
    num_events_generated = 0
    while num_TT < num_TT_events_to_generate:
        # convention for this test:
        # - l0=el, l1=mu
        # - l0_pt > l1_pt
        num_events_generated += 1
        l0_is_real = random_is_real(fraction_of_real_l0)
        l1_is_real = random_is_real(fraction_of_real_l1)
        l0_is_fake = not l0_is_real
        l1_is_fake = not l1_is_real
        l0_npt = random()
        l1_npt = random()
        l0_npt, l1_npt = (l0_npt, l1_npt) if l0_npt>l1_npt else (l1_npt, l0_npt)
        l0_pt  = pt_min + l0_npt*(pt_max - pt_min)
        l1_pt  = pt_min + l1_npt*(pt_max - pt_min)
        l0_ipt = int(l0_npt*n_pt_bins)
        l1_ipt = int(l1_npt*n_pt_bins)
        l0_is_tight = random_is_tight(eff_el_real[l0_ipt] if l0_is_real else eff_el_fake[l0_ipt])
        l1_is_tight = random_is_tight(eff_mu_real[l1_ipt] if l1_is_real else eff_mu_fake[l1_ipt])
        l0_is_loose = not l0_is_tight
        l1_is_loose = not l1_is_tight
        RR = (l0_is_real and l1_is_real)
        RF = (l0_is_real and l1_is_fake)
        FR = (l0_is_fake and l1_is_real)
        FF = (l0_is_fake and l1_is_fake)
        TT = (l0_is_tight and l1_is_tight)
        TL = (l0_is_tight and l1_is_loose)
        LT = (l0_is_loose and l1_is_tight)
        LL = (l0_is_loose and l1_is_loose)
        num_RR += 1 if RR else 0
        num_RF += 1 if RF else 0
        num_FR += 1 if FR else 0
        num_FF += 1 if FF else 0
        num_TT += 1 if TT else 0
        num_TL += 1 if TL else 0
        num_LT += 1 if LT else 0
        num_LL += 1 if LL else 0
        num_FF_TT += 1 if (FF and TT) else 0
        num_FR_TT += 1 if (FR and TT) else 0
        num_RF_TT += 1 if (RF and TT) else 0
        num_RR_TT += 1 if (RR and TT) else 0
        num_FF_LL += 1 if FF else 0
        num_FR_LL += 1 if FR else 0
        num_RF_LL += 1 if RF else 0
        num_RR_LL += 1 if RR else 0
        # if l0_is_real : m_r0.count(0, l0_is_tight)
        # else          : m_f0.count(0, l0_is_tight)
        # if l1_is_real : m_r1.count(0, l1_is_tight)
        # else          : m_f1.count(0, l1_is_tight)

        l1 = l1.isEl(True ).isTight(l0_is_tight).pt(l0_pt)#.eta(0.1)
        l0 = l0.isEl(False).isTight(l1_is_tight).pt(l1_pt)#.eta(0.1)
        weight = matrix.getTotalFake(l0, l1, region_index, systematic)
        tot_fake_weight += weight

    print "Summary of the test:"
    print "generated: {0} events".format(num_events_generated)
    print "of which {0} were TT".format(num_TT)
    print "of which {0} were due to fake ({1} FF + {2} FR + {3} RF)".format(num_FF_TT + num_FR_TT + num_RF_TT,
                                                                            num_FF_TT, num_FR_TT, num_RF_TT)
    print "estimated number of fakes: {0:.3f}".format(tot_fake_weight)
    print "\n"

def matrix_from_efficiencies(r0, r1, f0, f1):
    "using the same convention as in ATL-COM-PHYS-2011-144"
    m = r.TMatrixD(4, 4)
    m[0][0], m[0][1], m[0][2], m[0][3] = r0*r1,           r0*f1,           f0*r1,           f0*f1
    m[1][0], m[1][1], m[1][2], m[1][3] = r0*(1.-r1),      r0*(1.-f1),      f0*(1.-r1),      f0*(1.-f1)
    m[2][0], m[2][1], m[2][2], m[2][3] = (1.-r0)*r1,      (1.-r0)*f1,      (1.-f0)*r1,      (1.-f0)*f1
    m[3][0], m[3][1], m[3][2], m[3][3] = (1.-r0)*(1.-r1), (1.-r0)*(1.-f1), (1.-f0)*(1.-r1), (1.-f0)*(1.-f1)
    return m

def matrix_selected_events(n_tt, n_tl, n_lt, n_ll):
    "using the same convention as in ATL-COM-PHYS-2011-144"
    m = r.TMatrixD(4, 1)
    m[0][0] = n_tt
    m[1][0] = n_tl
    m[2][0] = n_lt
    m[3][0] = n_ll
    return m

def matrix_LL_events_compositions(n_RR_LL, n_RF_LL, n_FR_LL, n_FF_LL):
    "using the same convention as in ATL-COM-PHYS-2011-144"
    m = r.TMatrixD(4, 1)
    m[0][0] = n_RR_LL
    m[1][0] = n_RF_LL
    m[2][0] = n_FR_LL
    m[3][0] = n_FF_LL
    return m

def n_fake_from_efficiencies_and_selected(r0, r1, f0, f1, m_sel):
    "number of TT events due to fake"
    m_eff = matrix_from_efficiencies(r0, r1, f0, f1)
    m_inverted_eff = r.TMatrixD(m_eff)
    m_inverted_eff.Invert()
    m_LL = m_inverted_eff * m_sel
    n_RR_LL = m_LL[0][0]
    n_RF_LL = m_LL[1][0]
    n_FR_LL = m_LL[2][0]
    n_FF_LL = m_LL[3][0]
    n_RF_TT = r0*f1*n_RF_LL
    n_FR_TT = f0*r1*n_FR_LL
    n_FF_TT = f0*f1*n_FF_LL
    return n_RF_TT + n_FR_TT + n_FF_TT

#___________________________________________________________

main()
