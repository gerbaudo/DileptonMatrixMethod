#!/bin/env python

# assuming pt=x, eta=y

import datetime
import os
import sys
import ROOT as r
r.gROOT.SetBatch(1)

def main():
    input_filename = sys.argv[1]
    input_file = r.TFile.Open(input_filename)
    # input_file.ls()
    histo2d_names = ['el_real_eff2d_emuInc',
                     'mu_real_eff2d_emuInc',
                     'mu_fake_rate2d_emuInc',
                     'el_fake_rate2d_emuInc']
    titles = dict(zip(histo2d_names, ['Efficiency real electron', 'Efficiency real muon', 'Efficiency fake muon', 'Efficiency fake electron']))
    for hn in histo2d_names:
        h = input_file.Get(hn)
        histos1d = slice_histo(h)
        graphs = shifted_graphs(histos1d)
        pad_master = histos1d[0].Clone()
        pad_master.Reset()
        title = "Canvas created on {} from {} using {}".format(datetime.date.today().isoformat(),
                                                               os.getcwd(),
                                                               ' '.join(sys.argv))
        c = r.TCanvas(hn, title) # use c's title to encode metadata in the output  eps
        c.cd()
        c.SetRightMargin(1.5*c.GetRightMargin())
        pad_master.SetStats(False)
        xax, yax = pad_master.GetXaxis(), pad_master.GetYaxis()
        xax.SetLabelSize(yax.GetLabelSize())
        xax.SetTitleSize(yax.GetTitleSize())
        xax.SetTitle('p_{T} [GeV]')
        yax.SetTitle('#varepsilon(T|L)')
        pad_master.SetTitle(titles[hn])
        pad_master.Draw()
        for g in graphs:
            g.Draw('ep')
        draw_vertical_separators(pad_master, c)
        labels = eta_labels(h)
        leg = r.TLegend(0.8, 0.8, 1.0, 1.0)
        leg.SetFillColor(0)
        colors = [r.kRed, r.kBlue, r.kBlack]
        markers = [r.kStar, r.kCircle, r.kMultiply]
        for gr, lab, color, marker in zip(graphs, labels, colors, markers):
            gr.SetLineColor(color)
            gr.SetMarkerColor(gr.GetLineColor())
            gr.SetMarkerStyle(marker)
            leg.AddEntry(gr, lab, 'p')
        leg.Draw()
        c._leg = leg
        c.Update()
        c.SaveAs("/tmp/{}.png".format(hn))
        c.SaveAs("/tmp/{}.eps".format(hn))
        
def slice_histo(h2d):
    histos = []
    for i in range(1, 1+h2d.GetNbinsY()):
        histos = [h2d.ProjectionX(h2d.GetName()+"_bin%d"%i, i, i, 'e')
                  for i in range(1, 1+h2d.GetNbinsY())]
    return histos

def shifted_graphs(histos):
    "given N histos with overlapping markers, build N graphs with a width/N offset from one to another"
    assert len(set([h.GetNbinsX() for h in histos]))==1,"shifted_graphs: n bins ={}, expected {} values".format(set(h.GetNbinsX() for h in histos), 1)
    graphs = []
    for ih, h in enumerate(histos):
        gr = r.TGraphAsymmErrors(h)
        gr.SetName(h.GetName()+'_gr')
        for ib in range(1, 1+h.GetNbinsX()):
            bin_width = h.GetBinWidth(ib)
            xlo = h.GetBinLowEdge(ib)
            xhi = xlo+bin_width
            xcenter = 0.5*(xlo+xhi)
            step = bin_width/float(len(histos))
            offset = 0.5*step
            xcenter = xlo+offset+(ih)*step
            ycenter = h.GetBinContent(ib)
            erry = h.GetBinError(ib)
            gr.SetPoint(ib-1, xcenter, ycenter)
            errxlo, errxhi = xcenter-xlo, xhi-xcenter
            errxlo, errxhi = 0.0, 0.0
            gr.SetPointError(ib-1, errxlo, errxhi, erry, erry)
        graphs.append(gr)
    return graphs

def referenceLine(xmin=0., xmax=100.0, ymin=1.0, ymax=1.0) :
    l1 = r.TLine(xmin, ymin, xmax, ymax)
    l1.SetLineStyle(3)
    l1.SetLineColor(r.kGray+1)
    return l1

def draw_vertical_separators(h1d, canvas):
    ymin, ymax = 0.0, 1.05 #h1d.GetMinimum(), h1d.GetMaximum()
    lines = []
    for ib in range(1+1, 1+h1d.GetNbinsX()):
        x = h1d.GetBinLowEdge(ib)
        l = referenceLine(x, x, ymin, ymax)
        l.Draw()
        lines.append(l)
    canvas._lines = lines

def eta_labels(h2d):
    eta_axis = h2d.GetYaxis()
    return ["%.2f<|#eta|<%.2f"%(eta_axis.GetBinLowEdge(b), eta_axis.GetBinLowEdge(b)+eta_axis.GetBinWidth(b))
            for b in range(1, 1+eta_axis.GetNbins())]
if __name__=='__main__':
    main()
