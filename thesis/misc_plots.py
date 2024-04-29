# make plots of distribution and projections. 
# make plots comparing fits next to each other. make plot of p_T and 2D distributions for section on histogramming
# Notice a cut in eta invariably leads to a cut in delta eta; however a cut in delta eta is not the same as a cut in eta. 

import ROOT as rt
import numpy as np
import matplotlib.pyplot as plt

rt.gStyle.SetOptFit(1111) # show fit panel w/ fit prob, chi2/ndf, fit params, errors

import sys
sys.path.append('../classes/')
sys.path.append('../utils/')

from SparseAnalyzer import SparseAnalyzer #, FitFunction, FitRange, CurveFit
from formatting import process_pave, process_tgraph, process_hist

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
    rt.gStyle.SetOptStat(11)
    canvas = rt.TCanvas(f'{sparse.name}_plots_canvas', f'{sparse.name}_plots')
    canvas.SetCanvasSize(1800, 800)
    canvas.Divide(3, 2)

    canvas.cd(1)
    process_hist(trigEta)
    trigEta.GetYaxis().SetRangeUser(0, 1.3*trigEta.GetMaximum())
    trigEta.SetTitle('Trigger #eta')
    trigEta.SetXTitle('#eta')
    trigEta.SetYTitle('counts')
    trigEta.Draw('COLZ')


    canvas.cd(2)
    process_hist(assocEta)
    assocEta.SetTitle('Associate #eta')
    assocEta.SetXTitle('#eta')
    assocEta.SetYTitle('counts')
    assocEta.Draw('COLZ')


    canvas.cd(3)
    process_hist(trigPt)
    trigPt.SetTitle('Trigger p_{T}')
    trigPt.SetXTitle('p_{T}')
    trigPt.SetYTitle('counts')
    trigPt.Draw('COLZ')


    canvas.cd(4)
    process_hist(assocPt)
    assocPt.SetTitle('Associate p_{T}')
    assocPt.SetXTitle('p_{T}')
    assocPt.SetYTitle('counts')
    assocPt.Draw('COLZ')


    canvas.cd(5)
    process_hist(dPhiDEta)
    dPhiDEta.SetTitle('Correlation')
    dPhiDEta.SetXTitle('#Delta#varphi')
    dPhiDEta.SetYTitle('#Delta#eta')
    dPhiDEta.GetYaxis().SetRangeUser(-1.2, 1.199)
    dPhiDEta.Draw('SURF1')


    canvas.cd(6)
    process_hist(dEta)
    dEta.SetTitle('#Delta#eta')
    dEta.SetXTitle('#Delta#eta')
    dEta.SetYTitle('counts')
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

charged_hadron_pt = chargedHadronDist.Project3D('x')
charged_hadron_pt.Rebin(2)
charged_hadron_pt.GetYaxis().SetRangeUser(0, 1.3 * charged_hadron_pt.GetMaximum())
charged_hadron_pt.SetName('charged_hadron_pt')

pt_canvas.Draw()
charged_hadron_pt.Draw()
pt_canvas.SaveAs('charged_hadron_pt.pdf')

#########################
# 1D and 2D distributions
#########################

pt_eta_dist = chargedHadronDist.Project3D('xy')
pt_phi_dist = chargedHadronDist.Project3D('xz')

two_d_canvas = rt.TCanvas()
two_d_canvas.Divide(1, 2)

two_d_canvas.cd(1)
pt_eta_dist.Draw('SURF1')

two_d_canvas.cd(2)
pt_phi_dist.Draw('SURF1')

two_d_canvas.SaveAs('two_d_dists.pdf')


