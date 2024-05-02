# make plots of distribution and projections. 
# make plots comparing fits next to each other. make plot of p_T and 2D distributions for section on histogramming
# Notice a cut in eta invariably leads to a cut in delta eta; however a cut in delta eta is not the same as a cut in eta. 

import ROOT as rt
import numpy as np
import matplotlib.pyplot as plt

rt.gStyle.SetOptFit(1111) # show fit panel w/ fit prob, chi2/ndf, fit params, errors
rt.gStyle.SetOptStat(0)

import sys
sys.path.append('../classes/')
sys.path.append('../utils/')

from SparseAnalyzer import SparseAnalyzer #, FitFunction, FitRange, CurveFit
from formatting import process_pave, process_tgraph, process_hist, process_hist2D

from CurveFit import FitFunction, FitRange, CurveFit

EPSILON = 0.0001

# load the ROOT file
inFile = rt.TFile.Open("../dists/analysisMultiFine.root")

# load single particle histograms, TH3D
#   axes are: pT, phi, eta
allParticleDist = inFile.allParticleDist
chargedHadronDist = inFile.chargedHadronDist
lambdaDist = inFile.lambdaDist
triggeredLambdaDist = inFile.triggeredLambdaDist

# load correlation dists, 6D THnSparse 
#   axes are: trigger eta, associated eta, trigger pt, associated pt, delta phi, delta eta
hlDist = inFile.hLambdaDist
hhDist = inFile.hhDist


hh = SparseAnalyzer(hhDist, name='hh', apply_cut=False)
hh08 = SparseAnalyzer(hhDist, name='hh08', etaAssoc=[-0.8, 0.8-EPSILON], apply_cut=True)
hh12 = SparseAnalyzer(hhDist, name='hh12', etaAssoc=[-1.2, 1.2-EPSILON], apply_cut=True)
hh20 = SparseAnalyzer(hhDist, name='hh20', etaAssoc=[-2.0, 2.0-EPSILON], apply_cut=True)

hl = SparseAnalyzer(hlDist, name='hl', apply_cut=False)
hl08 = SparseAnalyzer(hlDist, name='hl08', etaAssoc=[-0.8, 0.8-EPSILON], apply_cut=True)
hl12 = SparseAnalyzer(hlDist, name='hl12', etaAssoc=[-1.2, 1.2-EPSILON], apply_cut=True)
hl20 = SparseAnalyzer(hlDist, name='hl20', etaAssoc=[-2.0, 2.0-EPSILON], apply_cut=True)

distributions = [hh, hh08, hh12, hh20, hl, hl08, hl12, hl20]

def make_projection_plots(sparse):
    '''
    Makes 1D histograms of quantities from THnSparse, and 2D delta phi delta eta correlation. 
    No single particle cuts applied in this function. 

    Args
        sparse (SparseAnalyzer):
    
    Returns
        canvas (TCanvas): 
    '''

    # project to make 1D hists
    trigEta = sparse.make_1D_hist('trigger eta')
    assocEta = sparse.make_1D_hist('associate eta')
    trigPt = sparse.make_1D_hist('trigger pT')
    assocPt = sparse.make_1D_hist('associate pT')
    dEta = sparse.make_1D_hist('delta eta') 

    # project into 2D correlations and delta eta dist
    dPhiDEta = sparse.make_2D_correlation()
    #dEta = dPhiDEta.ProjectionY()

    # make canvas
    # rt.gStyle.SetOptStat(11)
    canvas = rt.TCanvas(f'{sparse.name}_plots_canvas', f'{sparse.name}_plots')
    canvas.SetCanvasSize(1800, 800)
    canvas.Divide(3, 2)

    canvas.cd(1)
    rt.gPad.SetLeftMargin(0.2)
    rt.gPad.SetBottomMargin(0.16)
    process_hist(trigEta)
    trigEta.GetYaxis().SetRangeUser(0, 1.3*trigEta.GetMaximum())
    trigEta.SetTitle('Trigger #eta Distribution')
    trigEta.SetXTitle('#eta')
    trigEta.SetYTitle('counts')
    trigEta.GetYaxis().SetTitleOffset(1.5)
    trigEta.Draw('COLZ')


    canvas.cd(2)
    rt.gPad.SetLeftMargin(0.2)
    rt.gPad.SetBottomMargin(0.16)
    process_hist(assocEta)
    assocEta.SetTitle('Associated #eta Distribution')
    assocEta.SetXTitle('#eta')
    assocEta.SetYTitle('counts')
    assocEta.GetYaxis().SetTitleOffset(1.5)
    assocEta.Draw('COLZ')


    canvas.cd(3)
    rt.gPad.SetLeftMargin(0.2)
    rt.gPad.SetBottomMargin(0.16)
    process_hist(trigPt)
    trigPt.SetTitle('Trigger p_{T} Distribution')
    trigPt.SetXTitle('p_{T} (GeV)')
    trigPt.SetYTitle('counts')
    trigPt.GetYaxis().SetTitleOffset(1.5)
    trigPt.Draw('COLZ')


    canvas.cd(4)
    rt.gPad.SetLeftMargin(0.2)
    rt.gPad.SetBottomMargin(0.16)
    process_hist(assocPt)
    assocPt.SetTitle('Associated p_{T} Distribution')
    assocPt.SetXTitle('p_{T} (GeV)')
    assocPt.SetYTitle('counts')
    assocPt.GetYaxis().SetTitleOffset(1.5)
    assocPt.Draw('COLZ')


    canvas.cd(5)
    rt.gPad.SetLeftMargin(0.2)
    rt.gPad.SetBottomMargin(0.16)
    process_hist(dPhiDEta)
    dPhiDEta.SetTitle('2D Angular Correlation')
    dPhiDEta.SetXTitle('#Delta#varphi (rad)')
    dPhiDEta.SetYTitle('#Delta#eta')
    dPhiDEta.GetYaxis().SetRangeUser(-1.2, 1.199)
    dPhiDEta.GetYaxis().SetTitleOffset(1.5)
    dPhiDEta.Draw('SURF1')


    canvas.cd(6)
    rt.gPad.SetLeftMargin(0.2)
    rt.gPad.SetBottomMargin(0.16)
    process_hist(dEta)
    dEta.SetTitle('#Delta#eta Distribution')
    dEta.SetXTitle('#Delta#eta')
    dEta.SetYTitle('counts')
    dEta.GetYaxis().SetTitleOffset(1.5)
    dEta.Draw()    

    return canvas

##################
# Projection Plots
##################

hh08_prjxn_canvas = make_projection_plots(hh08)
hl08_prjxn_canvas = make_projection_plots(hl08)

hh08_prjxn_canvas.SaveAs('hh08_projections.pdf')
hl08_prjxn_canvas.SaveAs('hl08_projections.pdf')


################
# pT, no eta cut
################

pt_canvas = rt.TCanvas()
pt_canvas.SetLeftMargin(0.14)

charged_hadron_pt = chargedHadronDist.Project3D('x')
charged_hadron_pt.Rebin(2)
charged_hadron_pt.GetYaxis().SetRangeUser(0, 1.3 * charged_hadron_pt.GetMaximum())
charged_hadron_pt.GetXaxis().SetRangeUser(0, 6)

charged_hadron_pt.SetTitle('Charged Hadron p_{T} (no acceptance cut)')
charged_hadron_pt.GetXaxis().SetTitle('p_{T} (GeV/c)')
charged_hadron_pt.GetYaxis().SetTitle('#frac{dN}{dp_{T}} (1/GeV)')
charged_hadron_pt.GetYaxis().SetTitleOffset(1.8)

pt_canvas.Draw()
charged_hadron_pt.Draw()
pt_canvas.SaveAs('charged_hadron_pt.pdf')

##################
# eta distribution
##################

eta_canvas = rt.TCanvas()
eta_canvas.SetLeftMargin(0.14)

charged_hadron_eta = chargedHadronDist.Project3D('z')
#charged_hadron_eta.Rebin(2)
charged_hadron_eta.GetYaxis().SetRangeUser(0, 1.3 * charged_hadron_eta.GetMaximum())
charged_hadron_eta.GetXaxis().SetRangeUser(-6, 6)

charged_hadron_eta.SetTitle('Charged Hadron #eta (no acceptance cut)')
charged_hadron_eta.GetXaxis().SetTitle('#eta')
charged_hadron_eta.GetYaxis().SetTitle('#frac{dN}{d#eta}')
charged_hadron_eta.GetYaxis().SetTitleOffset(1.8)

eta_canvas.Draw()
charged_hadron_eta.Draw()
eta_canvas.SaveAs('charged_hadron_eta.pdf')

##################
# phi distribution
##################

phi_canvas = rt.TCanvas()
phi_canvas.SetLeftMargin(0.14)

charged_hadron_phi = chargedHadronDist.Project3D('y')
#charged_hadron_phi.Rebin(2)
charged_hadron_phi.GetYaxis().SetRangeUser(0, 1.3 * charged_hadron_phi.GetMaximum())
charged_hadron_phi.GetXaxis().SetRangeUser(0, 2 * rt.TMath.Pi() - 0.01)

charged_hadron_phi.SetTitle('Charged Hadron #varphi (no acceptance cut)')
charged_hadron_phi.GetXaxis().SetTitle('#varphi (rad)')
charged_hadron_phi.GetYaxis().SetTitle('#frac{dN}{d#varphi} (1/rad)')
charged_hadron_phi.GetYaxis().SetTitleOffset(1.8)

phi_canvas.Draw()
charged_hadron_phi.Draw()
phi_canvas.SaveAs('charged_hadron_phi.pdf')


#########################
# 1D and 2D distributions
#########################

pt_eta_dist = chargedHadronDist.Project3D('yx')
pt_phi_dist = chargedHadronDist.Project3D('xz')

pt_eta_dist.Rebin2D(2, 2)
pt_phi_dist.Rebin2D(2, 2)

two_d_canvas = rt.TCanvas('twod', 'twod', 1400, 700)

two_d_canvas.Divide(2, 1)
title_pad = rt.TPad("all", "all", 0, 0, 1, 1)
title_pad.SetFillStyle(4000) # transparent
title_pad.Draw()

rt.gStyle.SetOptTitle(0)

two_d_canvas.cd(1)
two_d_canvas.cd(1).SetPhi(225.)
rt.gPad.SetLeftMargin(0.2)
rt.gPad.SetBottomMargin(0.16)
#two_d_canvas.cd(1).SetLeftMargin(0.)

pt_eta_dist.SetXTitle('p_{T} (GeV)')
pt_eta_dist.SetYTitle('#varphi (rad)')
pt_eta_dist.SetZTitle('#frac{dN}{dp_{T} d#eta} (1/GeV)')
pt_eta_dist.GetZaxis().SetTitleOffset(2.2)
pt_eta_dist.GetXaxis().SetTitleOffset(2.2)
pt_eta_dist.GetYaxis().SetTitleOffset(2.2)
process_hist2D(pt_eta_dist)
pt_eta_dist.Draw('SURF1')

two_d_canvas.cd(2)
two_d_canvas.cd(2).SetPhi(225.)
rt.gPad.SetLeftMargin(0.2)
rt.gPad.SetBottomMargin(0.16)

pt_phi_dist.SetXTitle('#eta')
pt_phi_dist.SetYTitle('p_{T} (GeV)')
pt_phi_dist.SetZTitle('#frac{dN}{dp_{T}d#eta} (1/GeV)')

pt_phi_dist.GetZaxis().SetTitleOffset(2.2)

process_hist2D(pt_phi_dist)
pt_phi_dist.Draw('SURF1')

title_pad.cd()
latex = rt.TLatex()
latex.DrawLatexNDC(0.25, 0.92, "Charged Hadron Distributions (no acceptance cut)")

two_d_canvas.SaveAs('two_d_dists.pdf')


