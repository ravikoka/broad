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

widthRatiosPave = rt.TPaveText(0.5, 0.760593, 0.899425, 0.860169, 'NDC')
widthRatiosPave.AddText('pp, PYTHIA6, #sqrt{s}=14 TeV')
widthRatiosPave.AddText('4 < p_{T}^{trig} < 8 GeV/c, 2 < p_{T}^{assoc} < 4 GeV/c')
widthRatiosPave.AddText('|#eta^{trig}| < 0.8')
process_pave(widthRatiosPave, size=0.06)

def make_away_side_ratio_plots(sparse1, sparse2, eta_assoc_ranges):
    '''
    Calculate and graph away-side ratios for a range of associate eta cuts. 

    Args
        sparse1 (THnSparse): numerator
        sparse2 (THnSparse): denominator 

    Returns
        graph (TGraphErrors)
    '''
    
    yieldRatios = np.zeros(len(eta_assoc_ranges), dtype='d')
    yieldRatioUncerts = np.zeros(len(eta_assoc_ranges), dtype='d')
    etaLabels = np.array([np.abs(etaRange[0]) for etaRange in eta_assoc_ranges], dtype='d')

    for i, eta_assoc in enumerate(eta_assoc_ranges): 
        sparse1Cut = SparseAnalyzer(sparse1, 'numerator_dist', etaAssoc=eta_assoc)
        sparse2Cut = SparseAnalyzer(sparse2, 'denom_dist', etaAssoc=eta_assoc)

        numerator = sparse1Cut.integrate_away_side()
        denominator = sparse2Cut.integrate_away_side()

        yieldRatio = numerator / denominator
        yieldRatios[i] = yieldRatio

        sigma = yieldRatio * np.sqrt(1 / numerator + 1 / denominator)
        yieldRatioUncerts[i] = sigma

    etaLabels.ravel()
    etaLabelsUncerts = np.zeros(len(etaLabels), dtype='d').ravel()
    yieldRatios.ravel()
    yieldRatioUncerts.ravel()

    graph = rt.TGraphErrors(len(etaLabels), etaLabels, yieldRatios, ex=etaLabelsUncerts, ey=yieldRatioUncerts)

    return graph


ratioCanvas = rt.TCanvas()
rt.gPad.SetLeftMargin(0.19)
rt.gPad.SetBottomMargin(0.13)
#ratioCanvas.SetCanvasSize(800, 400)

ratioGraph = make_away_side_ratio_plots(hlDist, hhDist, eta_assoc_ranges=[[-2, 2-EPSILON], [-1.2, 1.2-EPSILON], [-0.8, 0.8-EPSILON]])
ratioGraph.SetMarkerStyle(107)
ratioGraph.SetTitle('Away Side Yield Ratios per #eta Cut')
ratioGraph.GetXaxis().SetTitle('|#eta^{assoc}| < x')
ratioGraph.GetYaxis().SetTitle('Y^{h-#Lambda}_{AS} / Y^{h-h}_{AS}')

process_tgraph(ratioGraph)
ratioGraph.GetYaxis().SetTitleOffset(1.0)

ratioGraph.SetMaximum(22e-3)
ratioGraph.SetMinimum(16e-3)

# ratioGraph.GetYaxis().SetRangeUser(18, 21)

ratioGraph.Draw('ap')

widthRatiosPave.Draw()
ratioCanvas.Draw()
ratioCanvas.SaveAs('away_yield_ratios.pdf')


denseRatioCanvas = rt.TCanvas()
denseRatioCanvas.SetCanvasSize(800, 400)

dense_eta = np.array([[-i*10**-1 , i*10**-1 - EPSILON] for i in range(8, 22, 2)])
denseRatioGraph = make_away_side_ratio_plots(hlDist, hhDist, eta_assoc_ranges=dense_eta)

linear = rt.TF1('linear_fit', 'pol1(0)')
denseRatioGraph.Fit('linear_fit')

denseRatioGraph.SetMarkerStyle(107)
denseRatioGraph.SetTitle('#(h - \Lambda)_{away} / (h - h)_{away}; |#eta|; away side ratio')

denseRatioGraph.Draw('ap')
denseRatioCanvas.Draw()
denseRatioCanvas.SaveAs('dense_yield_away.pdf')


def make_near_side_ratio_plots(sparse1, sparse2, eta_assoc_ranges):
    '''
    Calculate and graph near-side ratios for a range of associate eta cuts. 

    Args
        sparse1 (THnSparse): numerator
        sparse2 (THnSparse): denominator 

    Returns
        graph (TGraphErrors)
    '''
    
    yieldRatios = np.zeros(len(eta_assoc_ranges), dtype='d')
    yieldRatioUncerts = np.zeros(len(eta_assoc_ranges), dtype='d')
    etaLabels = np.array([np.abs(etaRange[0]) for etaRange in eta_assoc_ranges], dtype='d')

    for i, eta_assoc in enumerate(eta_assoc_ranges): 
        sparse1Cut = SparseAnalyzer(sparse1, 'numerator_dist', etaAssoc=eta_assoc)
        sparse2Cut = SparseAnalyzer(sparse2, 'denom_dist', etaAssoc=eta_assoc)

        numerator = sparse1Cut.integrate_near_side()
        denominator = sparse2Cut.integrate_near_side()

        yieldRatio = numerator / denominator
        yieldRatios[i] = yieldRatio

        sigma = yieldRatio * np.sqrt(1 / numerator + 1 / denominator)
        yieldRatioUncerts[i] = sigma

    etaLabels.ravel()
    etaLabelsUncerts = np.zeros(len(etaLabels), dtype='d').ravel()
    yieldRatios.ravel()
    yieldRatioUncerts.ravel()

    graph = rt.TGraphErrors(len(etaLabels), etaLabels, yieldRatios, ex=etaLabelsUncerts, ey=yieldRatioUncerts)

    return graph


nsYieldRatioCanvas = rt.TCanvas()
rt.gPad.SetLeftMargin(0.19)
rt.gPad.SetBottomMargin(0.13)

nsYieldRatioGraph = make_near_side_ratio_plots(hlDist, hhDist, eta_assoc_ranges=[[-2, 2-EPSILON], [-1.2, 1.2-EPSILON], [-0.8, 0.8-EPSILON]])
nsYieldRatioGraph.SetMarkerStyle(21)
nsYieldRatioGraph.SetMarkerColor(4)
nsYieldRatioGraph.SetLineColor(4)

nsYieldRatioGraph.SetTitle('Near Side Yield Ratios per #eta Cut')
nsYieldRatioGraph.GetXaxis().SetTitle('|#eta^{assoc}| < x')
nsYieldRatioGraph.GetYaxis().SetTitle('Y^{h-#Lambda}_{NS} / Y^{h-h}_{NS}')

nsYieldRatioGraph.SetMinimum(3e-3)
nsYieldRatioGraph.SetMaximum(6e-3)

process_tgraph(nsYieldRatioGraph)

#nsYieldRatioGraph.GetYaxis().SetRangeUser(4, 6)

nsYieldRatioGraph.Draw('ap')
widthRatiosPave.Draw()
nsYieldRatioCanvas.Draw()
nsYieldRatioCanvas.SaveAs('near_yield_ratios.pdf')

nsDenseRatioCanvas = rt.TCanvas()
nsDenseRatioCanvas.SetCanvasSize(800, 400)

dense_eta = np.array([[-i*10**-1 , i*10**-1 - EPSILON] for i in range(8, 22, 2)])
nsDenseRatioGraph = make_near_side_ratio_plots(hlDist, hhDist, eta_assoc_ranges=dense_eta)

linear = rt.TF1('linear_fit_ns', 'pol1(0)')
nsDenseRatioGraph.Fit('linear_fit_ns')

nsDenseRatioGraph.SetMarkerStyle(107)
nsDenseRatioGraph.SetTitle('#(h - \Lambda)_{away} / (h - h)_{away}; |#eta|; near side ratio')

nsDenseRatioGraph.Draw('ap')
nsDenseRatioCanvas.Draw()
nsDenseRatioCanvas.SaveAs('dense_yield_near.pdf')



def make_away_over_near_plots(sparse, eta_assoc_ranges):
    
    yieldRatios = np.zeros(len(eta_assoc_ranges), dtype='d')
    yieldRatioUncerts = np.zeros(len(eta_assoc_ranges), dtype='d')
    etaLabels = np.array([np.abs(etaRange[0]) for etaRange in eta_assoc_ranges], dtype='d')

    for i, eta_assoc in enumerate(eta_assoc_ranges): 
        sparseCut = SparseAnalyzer(sparse, 'cut_sparse', etaAssoc=eta_assoc)

        numerator = sparseCut.integrate_away_side()
        denominator = sparseCut.integrate_near_side()

        yieldRatio = numerator / denominator
        yieldRatios[i] = yieldRatio

        sigma = yieldRatio * np.sqrt(1 / numerator + 1 / denominator)
        yieldRatioUncerts[i] = sigma

    etaLabels.ravel()
    etaLabelsUncerts = np.zeros(len(etaLabels), dtype='d').ravel()
    yieldRatios.ravel()
    yieldRatioUncerts.ravel()

    graph = rt.TGraphErrors(len(etaLabels), etaLabels, yieldRatios, ex=etaLabelsUncerts, ey=yieldRatioUncerts)

    return graph
    


hhRatios = make_away_over_near_plots(hhDist, eta_assoc_ranges=[[-2, 2-EPSILON], [-1.2, 1.2-EPSILON], [-0.8, 0.8-EPSILON]])
hhRatios.SetMarkerStyle(22)
hhRatios.SetTitle('h-h yield ratio per #eta cut')

hlRatios = make_away_over_near_plots(hlDist, eta_assoc_ranges=[[-2, 2-EPSILON], [-1.2, 1.2-EPSILON], [-0.8, 0.8-EPSILON]])
hlRatios.SetMarkerStyle(22)
hlRatios.SetMarkerColor(4)
hlRatios.SetTitle('h-#Lambda yield ratio per #eta cut')


sameRatioCanvas = rt.TCanvas()
rt.gPad.SetLeftMargin(0.19)
rt.gPad.SetBottomMargin(0.13)


sameRatioMG = rt.TMultiGraph('yieldmg', 'away-side yield ratios per #eta cut')

sameRatioMG.Add(hhRatios, 'p')
sameRatioMG.Add(hlRatios, 'p')

sameRatioMG.GetXaxis().SetTitle('|#eta^{assoc}| < x')
sameRatioMG.GetYaxis().SetTitle('Y^{AS} / Y^{NS}')

sameRatioMG.Draw()

sameRatioLegend = sameRatioCanvas.BuildLegend(0.232758, 0.705508, 0.53304, 0.91525)
sameRatioLegend.SetBorderSize(0)
sameRatioLegend.SetFillStyle(0)

sameRatioCanvas.Draw()
sameRatioCanvas.SaveAs('away-over-near_yields.pdf')