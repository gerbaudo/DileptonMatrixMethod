#!/bin/env python

# Given an input fake matrix, test it with some input events
#
# The fake matrix file is usually one from the 'data' dir.
# The input events can be provided either from stdin or txt file.
#
# davide.gerbaudo@gmail.com
# 2013-07-25

import optparse
import os
import sys

def main():
    options = parse_options()
    file_name = options.input
    matrix_name = options.matrix_file
    region = options.region
    use_param_pt = options.param_pt
    use_param_pt_eta = options.param_pt_eta
    verbose = options.verbose

    if options.print_example_input:
        print_example_input()
        return
    if verbose:
        print "Using the following options:"
        print "matrix : {0}".format(matrix_name)
        print "region : {0}".format(region)
        print "parametrization : {0}".format('pt' if use_param_pt else 'pt_eta')

    load_packages()
    matrix = r.susy.fake.DileptonMatrixMethod()
    regions = r.std.vector('std::string')()
    regions.push_back(region)
    parametrization = r.susy.fake.Parametrization.PT if use_param_pt else r.susy.fake.Parametrization.PT_ETA
    p = parametrization
    matrix.configure(matrix_name, regions, p, p, p, p)
    systematic = r.susy.fake.Systematic.SYS_NOM
    region_index = matrix.getIndexRegion(region)
    l0 = r.susy.fake.Lepton(False, False, 10.0, 1.0)
    l1 = r.susy.fake.Lepton(False, False, 10.0, 1.0)
    pt_unit = 1.0e3 if options.mev2gev else 1.0
    input_stream = open(file_name) if file_name else sys.stdin
    for line in input_stream.readlines():
        line = line.strip()
        if is_valid_line(line):
            if verbose: print line.strip()
            parameters = parse_input_line(line)
            l0El, l0Tight, l0pt, l0eta = parameters[ ::2]
            l1El, l1Tight, l1pt, l1eta = parameters[1::2]
            l0 = l0.isEl(l0El).isTight(l0Tight).pt(l0pt*pt_unit).eta(l0eta)
            l1 = l1.isEl(l1El).isTight(l1Tight).pt(l1pt*pt_unit).eta(l1eta)
            weight = matrix.getTotalFake(l0, l1, region_index, systematic)
            print "weight {0:10.6f} l0 {1} l1 {2}".format(weight, l0.str(), l1.str())
        elif verbose:
            print line



#___________________________________________________________
def default_region() : return 'emuInc'
def parse_options():
    usage = """usage: %prog [options]
    Example:
    %prog -i file.txt -r emuInc --param-pt -v
or
    %prog --print-example-input | %prog -r emuInc --param-pt-eta -v
    """
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-i', '--input', help='input file; default stdin')
    parser.add_option('-m', '--matrix-file', help='fake matrix file')
    parser.add_option('-r', '--region', help='signal region; use print_available_regions.py to determine the available ones')
    parser.add_option('--param-pt', action='store_true', default=False, help='use 1d parametrization (vs. pt). Default one.')
    parser.add_option('--param-pt-eta', action='store_true', default=False, help='use 2d parametrization (vs. pt and eta)')
    parser.add_option('--print-example-input', action='store_true', default=False, help='print example input to stdout')
    parser.add_option('--mev2gev', action='store_true', default=False,
                      help='assume the input matrix has pt in MeV (previous format, pt values are now in GeV)')
    parser.add_option('-v', '--verbose', action='store_true', default=False, help='print more details about what is going on')
    (options, args) = parser.parse_args()
    if options.input and not os.path.exists(options.input) : parser.error("invalid input '{0}'".format(options.input))
    if not options.matrix_file : options.matrix_file = default_input_matrix()
    if not os.path.exists(options.matrix_file) : parser.error("invalid matrix file '{0}'".format(options.matrix_file))
    if not options.region : options.region = default_region()
    if not options.param_pt and not options.param_pt_eta : options.param_pt = True
    return options
#___________________________________________________________
def rootcoredir():
    return os.environ['ROOTCOREDIR']

def import_root():
    import ROOT as r
    r.gROOT.SetBatch(1)
    r.PyConfig.IgnoreCommandLineOptions = True # don't let root steal your cmd-line options
    return r

r = import_root()
def load_packages():
    "load the DileptonMatrixMethod package, as RootCore's load_packages.C would do"
    r.gSystem.AddIncludePath('-I"'+rootcoredir()+'/include"')
    r.gSystem.Load('libDileptonMatrixMethod')

def default_input_matrix():
    return os.path.join(rootcoredir(), 'data/DileptonMatrixMethod/FakeMatrix_Jul_26.root')
#__________________________________________________________

def print_example_input():
    """
    Formatting: one event per line.
    Pt values are in GeV. Comment ('#') and empty lines are skipped.
    Each line should look like:
    <isElectron 1> <isElectron 2> <isTight 1> <isTight 2> <pt 1> <pt 2> <eta 1> <eta 2>
    """
    print """
# format:
# <isElectron 1> <isElectron 2> <isTight 1> <isTight 2> <pt 1> <pt 2> <eta 1> <eta 2>
# ee TT
1 1 1 1 35.0 25.0 1.5 0.9
# emu TL
1 0 1 0 36.0 24.0 1.4 1.0
# mumu  LL
0 0 0 0 37.0 23.0 1.3 1.1
"""

def is_valid_line(line):
    line = line.strip()
    num_expected_fields = 8
    return len(line)>0 and line[0]!='#' and len(line.split())==num_expected_fields

def parse_input_line(line):
    words = line.strip().split()
    isEl1, isEl2, isTight1, isTight2, pt1, pt2, eta1, eta2 = words
    isEl1, isEl2, isTight1, isTight2 = bool(int(isEl1)), bool(int(isEl2)), bool(int(isTight1)), bool(int(isTight2))
    pt1, pt2, eta1, eta2 = float(pt1), float(pt2), float(eta1), float(eta2)
    return isEl1, isEl2, isTight1, isTight2, pt1, pt2, eta1, eta2
#___________________________________________________________

if __name__=='__main__':
    main()
