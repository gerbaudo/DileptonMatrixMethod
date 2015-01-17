#!/bin/env python

# Add some reasonable syst uncertainties for the razor weighted
# average, based on what we had for the CR_SSInc1j weighted average.
# These uncertainties account for the stat uncertainty on the scale
# factor and for the uncertainty on the composition.

# davide.gerbaudo@gmail.com
# Jan 2015

import subprocess
from math import sqrt

import ROOT as r
r.gROOT.SetBatch(1)

from compare_matrices import binContents, getBinIndices, formatted_histo_contents

# avg values from mu_fake_rate2d_CR_SSInc1j
mu_fake_fractional_unc = [0.0, 0.5*(0.23+0.29), 0.5*(0.40+0.52)]
# avg values from el_fake_rate2d_CR_SSInc1j
el_fake_fractional_unc = [0.0, 0.5*(0.17+0.14), 0.5*(0.18+0.17)]

# avg values from  mu_fake_rate2d_CR_SSInc1j_frac_[up,do]
mu_comp_fractional_unc = [0.0, 0.25*(0.00+0.02+0.00+0.22), 0.25*(0.00+0.04+0.00+0.23)]
# avg values from  el_fake_rate2d_CR_SSInc1j_frac_[up,do]
el_comp_fractional_unc = [0.0, 0.25*(0.15+0.18+0.14+0.13), 0.25*(0.15+0.01+0.12+0.10)]

def main():
    input_fname_template = 'FinalFakeHist_May_20.root' # WH matrix file, from which we'll get the uncertainties
    input_fname = 'FakeMatrix_Nov_26.root'             # razor file, from which we'll get the nominal fake rates
    output_fname = 'FakeMatrix_Nov_26_with_syst.root'  # destination file
    # input_fname = 'FakeMatrix_Dec_03.root'             # razor file, from which we'll get the nominal fake rates
    # output_fname = 'FakeMatrix_Dec_03_with_syst.root'  # destination file

    print_existing_uncertainties(input_fname_template)

    clone_file(input_fname, output_fname)

    add_uncertainties(output_fname)

def print_frac_uncertainty(h):
    "given a 1D or a 2D histo, print its per-bin fractional uncertainty"
    print h.GetName(),' : bin errors'
    for b in getBinIndices(h):
        bc = h.GetBinContent(b)
        be = h.GetBinError(b)
        print "{}) : {:.3f} +/- {:.3f} ({:.0%})".format(b, bc, be, (be/bc) if bc else 0.0)
def print_frac_delta(h_n, h_u, h_d):
    "given two 1D or 2D histos, print its per-bin fractional diff"
    print h_n.GetName(),' : nominal, up, down'
    for b in getBinIndices(h_n):
        bcn = h_n.GetBinContent(b)
        bcu = h_u.GetBinContent(b)
        bcd = h_d.GetBinContent(b)
        print "{}) : {:.3f} + {:.3f} ({:.0%}) / - {:.3f} ({:.0%})".format(b, bcn,
                                                                          bcu,
                                                                          ((bcu-bcn)/bcn) if bcn else 0.0,
                                                                          bcd,
                                                                          ((bcd-bcn)/bcn) if bcn else 0.0)
def print_existing_uncertainties(filename):        
    input_file = r.TFile.Open(filename)
    input_file.cd()

    for hname in ['mu_fake_rate2d_CR_SSInc1j',
                  'el_fake_rate2d_CR_SSInc1j',
                  'el_fake_rate_CR_SSInc1j',
                  'mu_fake_rate_CR_SSInc1j']:
        h = input_file.Get(hname)
        print "{0} : {1}".format(hname, formatted_histo_contents(h))
        print_frac_uncertainty(h)

    mu_nom   = input_file.Get('mu_fake_rate2d_CR_SSInc1j')
    el_nom   = input_file.Get('el_fake_rate2d_CR_SSInc1j')
    mu_fr_up = input_file.Get('mu_fake_rate2d_CR_SSInc1j_frac_up')
    mu_fr_do = input_file.Get('mu_fake_rate2d_CR_SSInc1j_frac_do')
    el_fr_up = input_file.Get('el_fake_rate2d_CR_SSInc1j_frac_up')
    el_fr_do = input_file.Get('el_fake_rate2d_CR_SSInc1j_frac_do')

    print "composition uncertainties:"
    print '\n   muon_fake'
    print_frac_delta(mu_nom, mu_fr_up, mu_fr_do)
    
    print '\n   elec_fake'
    print_frac_delta(el_nom, el_fr_up, el_fr_do)
    input_file.Close()

def getCommandOutput(command):
    "lifted from supy (https://github.com/elaird/supy/blob/master/utils/io.py)"
    p = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    stdout,stderr = p.communicate()
    return {"stdout":stdout, "stderr":stderr, "returncode":p.returncode}


def clone_file(inputfilename, outputfilename):
    cmd = "cp -p {} {}".format(inputfilename, outputfilename)
    print cmd
    getCommandOutput(cmd)

def add_uncertainties(filename):
    input_file = r.TFile.Open(filename, 'update')
    input_file.cd()
    def set_scale_error(h, errors):
        print "set_scale_error for {}".format(h.GetName())
        nbins = h.GetNbinsX()
        assert nbins==len(errors),"invalid error array {}, expect {} elements".format(str(error), nbins)
        for i, err in enumerate(errors):
            bin = i+1
            bc = h.GetBinContent(bin)
            old_err = h.GetBinError(bin)
            new_contr = bc*err
            new_err = sqrt(old_err**2 + new_contr**2)
            print "[{}] error was {:.4f}, adding in quad {:.4f} -> {:.4f}".format(bin, old_err, new_contr, new_err)
            h.SetBinError(bin, new_err)
    el_nom   = input_file.Get('el_fake_rate_razor')
    mu_nom   = input_file.Get('mu_fake_rate_razor')
    set_scale_error(el_nom, el_fake_fractional_unc)
    set_scale_error(mu_nom, mu_fake_fractional_unc)


    mu_nom   = input_file.Get('mu_fake_rate_razor')
    el_nom   = input_file.Get('el_fake_rate_razor')
    mu_fr_up = mu_nom.Clone('mu_fake_rate_razor_frac_up')
    mu_fr_do = mu_nom.Clone('mu_fake_rate_razor_frac_do')
    el_fr_up = el_nom.Clone('el_fake_rate_razor_frac_up')
    el_fr_do = el_nom.Clone('el_fake_rate_razor_frac_do')

    def set_frac_variation(h, frac_variations, updown=+1.0):
        print "set_frac_variation for {}".format(h.GetName())
        nbins = h.GetNbinsX()
        assert nbins==len(frac_variations),"invalid variations array {}, expect {} elements".format(str(frac_variations), nbins)
        for i, frac_variation in enumerate(frac_variations):
            bin = i+1
            old_bc = h.GetBinContent(bin)
            new_bc = old_bc*(1.0+updown*frac_variation)
            print "[{}] rate was {:.4f}, now {:.4f} ({:.4f})".format(bin, old_bc, new_bc, updown*frac_variation)
            h.SetBinContent(bin, new_bc)

    set_frac_variation(mu_fr_up, mu_comp_fractional_unc, updown=+1.0)
    set_frac_variation(mu_fr_do, mu_comp_fractional_unc, updown=-1.0)
    set_frac_variation(el_fr_up, el_comp_fractional_unc, updown=+1.0)
    set_frac_variation(el_fr_do, el_comp_fractional_unc, updown=-1.0)
            
    input_file.Write()
    input_file.Close()
            

if __name__=='__main__':
    main()
