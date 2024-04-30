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

##########################
## Width Extraction
##########################

# standard gaussian fit
dPhiCanvas = rt.TCanvas('dphi_canvas', 'dphi_canvas')

dPhiCanvas.SetCanvasSize(2400, 800)
dPhiCanvas.Divide(4, 2)

guesses = [[15000, 0.2, 10000, 0.2], [20000, 1, 8000, 1], [20000, 1, 15000, 1], [20000, 1, 20000, 1],
           [200, 1, 800, 1], [100, 1, 150, 1], [100, 0.19, 250, 1], [70, 1, 400, 1]]
awayStDev = np.zeros(shape=(len(distributions), 2))
nearStDev = np.zeros(shape=(len(distributions), 2))
for i, (sparse, guess) in enumerate(zip(distributions, guesses)):
    dPhiCanvas.cd(i + 1)

    dPhi = sparse.make_delta_phi()
    dPhi, params = sparse.fit_delta_phi(dPhi, guess, ratio_plot=False)
    
    nearStDev[i] = params[2]
    awayStDev[i] = params[5]
    print(params[2], params[5])
    print(nearStDev)
    

    pave = sparse.make_pave()    
    dPhi.Draw()
    pave.Draw()
        
dPhiCanvas.SaveAs('delta_phi_plots.png')

#############################
## Width Ratios Calculations
#############################

# plot ratio of yields in near side and away side, for same type of correlation
widthsNearOverAwayCanvas = rt.TCanvas()
widthsNearOverAwayCanvas.SetCanvasSize(800, 400)

# get rid of ratios for distributions with no eta cuts
filter = [False, True, True, True, 
          False, True, True, True]

nearStDev = nearStDev[filter]
awayStDev = awayStDev[filter]

# get widths for h-h azimuthal correlations
hhAwayWidths = awayStDev[:3, 0]
hhNearWidths = nearStDev[:3, 0]
hhAwayErrs = awayStDev[:3, 1]
hhNearErrs = nearStDev[:3, 1]

# get widths for h-Lambda azimuthal correlations
hlAwayWidths = awayStDev[3:, 0]
hlNearWidths = nearStDev[3:, 0]
hlAwayErrs = awayStDev[3:, 1]
hlNearErrs = nearStDev[3:, 1]

# typical error propagation for ratio
errorProp = lambda x, xerr, y, yerr: (x / y) * np.sqrt((xerr / x)**2 + (yerr / y)**2)

# calc h-h near side width / h-h away side width
# hhRatios = hhNearWidths / hhAwayWidths
# hhErrs = errorProp(hhNearWidths, hhNearErrs, hhAwayWidths, hhAwayErrs)
hhRatios = hhAwayWidths / hhNearWidths
hhErrs = errorProp(hhAwayWidths, hhAwayErrs, hhNearWidths, hhNearErrs)

# calc h-Lambda near side width / h-h near side width
# hlRatios = hlNearWidths / hlAwayWidths
# hlErrs = errorProp(hlNearWidths, hlNearErrs, hlAwayWidths, hlAwayErrs)
hlRatios = hlAwayWidths / hlNearWidths
hlErrs = errorProp(hlAwayWidths, hlAwayErrs, hlNearWidths, hlNearErrs)

# calc h-Lambda away side width / h-h away side width
widthRatios = hlAwayWidths / hhAwayWidths
ratioErrs = errorProp(hlAwayWidths, hlAwayErrs, hhAwayWidths, hhAwayErrs)
# ratioErrs = widthRatios * np.sqrt(1 / hlAwayWidths + 1 / hhAwayWidths)

print(widthRatios, ratioErrs)

etaAxis = np.array([0.8, 1.2, 2.0], dtype='d')
etaErrs = np.zeros(len(etaAxis))


##########################
## Away-Side Width Ratios
##########################

widthRatiosCanvas = rt.TCanvas('widthRatiosCanvas', 'widthRatiosCanvas')
rt.gPad.SetLeftMargin(0.15)
rt.gPad.SetBottomMargin(0.15)

widthRatiosPave = rt.TPaveText(0.25, 0.65, 0.55, 0.85, 'NDC') 
widthRatiosPave.AddText('pp, PYTHIA6, #sqrt{s}=7 TeV')
widthRatiosPave.AddText('4 < p_{T}^{trig} < 8 GeV/c, 2 < p_{T}^{assoc} < 4 GeV/c')
widthRatiosPave.AddText('|#eta^{trig}| < 0.8')
process_pave(widthRatiosPave, size=0.04)#0.06)


widthRatiosGraph = rt.TGraphErrors(len(widthRatios), etaAxis, widthRatios, ex=etaErrs, ey=ratioErrs)
widthRatiosGraph.SetMarkerStyle(107)

widthRatiosGraph.SetTitle('Away-side width ratio for each #eta^{assoc} range')
#widthRatiosGraph.GetYaxis().SetTitle('Width Ratio #left(#frac{h-#Lambda}{h-h}#right)')
widthRatiosGraph.GetYaxis().SetTitle('#sigma^{h-#Lambda}_{AS} / #sigma^{h-h}_{AS}')
widthRatiosGraph.GetXaxis().SetTitle('|#eta^{assoc}| < x')

process_tgraph(widthRatiosGraph)

#widthRatiosGraph.GetYaxis().SetTitleOffset(1.1)
#widthRatiosGraph.GetYaxis().SetTitleSize(0.055)

widthRatiosGraph.GetYaxis().SetRangeUser(0.8, 1.4)#0.8, 1.2)
#widthRatiosGraph.GetYaxis().LabelsOption('v')

line = rt.TLine(0.68, 1.0, 2.12, 1.0)
line.SetLineColor(4)
line.SetLineStyle(2)
line.SetLineWidth(2)

widthRatiosGraph.Draw('ap')
widthRatiosPave.Draw()
line.Draw()
widthRatiosCanvas.Draw()
widthRatiosCanvas.SaveAs('away_side_width_ratios.pdf')


##########################
## Near-Side Width Ratios
##########################

nsWidthRatiosCanvas = rt.TCanvas('nsWidthRatiosCanvas', 'nsWidthRatiosCanvas')
rt.gPad.SetLeftMargin(0.15)
rt.gPad.SetBottomMargin(0.15)

nsWidthRatios = hlNearWidths / hhNearWidths
nsWidthErrs = errorProp(hlNearWidths, hlNearErrs, hhNearWidths, hhNearErrs)

nsWidthRatiosGraph = rt.TGraphErrors(len(etaAxis), etaAxis, nsWidthRatios, ex=etaErrs, ey=nsWidthErrs)
nsWidthRatiosGraph.SetMarkerStyle(21)
nsWidthRatiosGraph.SetMarkerColor(4)
nsWidthRatiosGraph.SetLineColor(4)

nsWidthRatiosGraph.SetTitle('Near-side Width Ratio for each #eta^{assoc} range')
nsWidthRatiosGraph.GetYaxis().SetTitle('#sigma^{h-#Lambda}_{NS} / #sigma^{h-h}_{NS}')
nsWidthRatiosGraph.GetXaxis().SetTitle('|#eta^{assoc}| < x')

process_tgraph(nsWidthRatiosGraph)
nsWidthRatiosGraph.SetMaximum(1.4)
nsWidthRatiosGraph.SetMinimum(0.8)

line.SetLineColor(rt.kBlack)

nsWidthRatiosGraph.Draw('ap')
widthRatiosPave.Draw()
line.Draw()
nsWidthRatiosCanvas.Draw()
nsWidthRatiosCanvas.SaveAs('near_side_width_ratios.pdf')


###############################
## Away-Over-Near Width Ratios
###############################

#hhRatiosCanvas = rt.TCanvas()
#rt.gPad.SetLeftMargin(0.19)
#rt.gPad.SetBottomMargin(0.13)

hhRatiosGraph = rt.TGraphErrors(len(hhRatios), etaAxis, hhRatios, ex=etaErrs, ey=hhErrs)
hhRatiosGraph.SetTitle('h-h Width Ratios per #eta Cut')
#hhRatiosGraph.GetYaxis().SetTitle('h-h Width Ratios #left(#frac{NS}{AS}#right)')
hhRatiosGraph.GetYaxis().SetTitle('#sigma^{h-h}_{AS} / #sigma^{h-h}_{NS}')
hhRatiosGraph.GetXaxis().SetTitle('|#eta^{assoc}| < x')

hhRatiosGraph.SetMarkerStyle(20)
hhRatiosGraph.SetMarkerColor(2)
hhRatiosGraph.SetLineColor(2)

#process_tgraph(hhRatiosGraph)

#hhRatiosGraph.SetMinimum(0)
#hhRatiosGraph.SetMaximum(2.3)

#hRatiosGraph.Draw('ap')
#hhRatiosPave.Draw()
#hhRatiosCanvas.Draw()

#hlRatiosCanvas = rt.TCanvas()
#rt.gPad.SetLeftMargin(0.19)
#rt.gPad.SetBottomMargin(0.13)

hlRatiosPave = widthRatiosPave.Clone()

hlRatiosGraph = rt.TGraphErrors(len(hlRatios), etaAxis, hlRatios, ex=etaErrs, ey=hlErrs)
hlRatiosGraph.SetTitle('h-#Lambda width ratios per #eta cut')
hlRatiosGraph.GetXaxis().SetTitle('|#eta^{assoc}| < x')
#hlRatiosGraph.GetYaxis().SetTitle('h-#Lambda Width Ratios #left(#frac{NS}{AS}#right)')
hlRatiosGraph.GetYaxis().SetTitle('#sigma^{h-#Lambda}_{AS} / #sigma^{h-#Lambda}_{NS}')

#process_tgraph(hlRatiosGraph)

hlRatiosGraph.SetMinimum(0.0)
hlRatiosGraph.SetMaximum(2.5)


hlRatiosGraph.SetMarkerStyle(21)
#hlRatiosGraph.Draw('ap')
#hlRatiosPave.Draw()
#hlRatiosCanvas.Draw()


compareRatiosCanvas = rt.TCanvas()
rt.gPad.SetLeftMargin(0.19)
rt.gPad.SetBottomMargin(0.13)

widthsMultiGraph = rt.TMultiGraph('widthmg', 'widthmg')

#hhRatiosGraph.SetMarkerColor(2)
#hlRatiosGraph.SetMarkerColor(1)

#hhRatiosGraph.SetMarkerStyle(8)
#hlRatiosGraph.SetMarkerStyle(21)


hhRatiosGraph.SetTitle('#sigma^{h-h}_{AS} / #sigma^{h-h}_{NS}')
hlRatiosGraph.SetTitle('#sigma^{h-#Lambda}_{AS} / #sigma^{h-#Lambda}_{NS}')
widthsMultiGraph.Add(hhRatiosGraph, 'p')
widthsMultiGraph.Add(hlRatiosGraph, 'p')


widthsMultiGraph.SetTitle('Comparison of Widths Ratios for each #eta^{assoc} range')
widthsMultiGraph.GetXaxis().SetTitle('|#eta^{assoc}| < x')
#widthsMultiGraph.GetYaxis().SetTitle('Width Ratios #left(#frac{NS}{AS}#right)')
widthsMultiGraph.GetYaxis().SetTitle('#sigma_{AS} / #sigma_{NS}')
widthsMultiGraph.GetYaxis().SetRangeUser(1.2, 2.6)#1.2, 2.5)

process_tgraph(widthsMultiGraph)

widthsMultiGraph.Draw('ap')

legend = compareRatiosCanvas.BuildLegend(0.571839, 0.68644, 0.87212, 0.89618)
legend.SetBorderSize(0)
legend.SetFillStyle(0)
legend.SetTextSize(0.04)#0.035)

widthsMultiGraphPave = rt.TPaveText(0.4, 0.3, 0.799425, 0.5, 'NDC')#0.189655, 0.764830, 0.589080, 0.86440677, 'NDC')
widthsMultiGraphPave.AddText('4 < p_{T}^{trig} < 8 GeV/c, 2 < p_{T}^{assoc} < 4 GeV/c')
widthsMultiGraphPave.AddText('|#eta^{trig}| < 0.8')
widthsMultiGraphPave = process_pave(widthRatiosPave, size=0.06)
widthsMultiGraphPave.SetTextSize(0.039)

# set coordinates of pave and legend, ROOT ignores the previous coords for reasons unknown to me
rt.gPad.Update()
widthsMultiGraphPave.SetX1NDC(0.29)#0.1896551724137931)
widthsMultiGraphPave.SetY1NDC(0.7)#0.760593220338983)
widthsMultiGraphPave.SetX2NDC(0.59)#0.5890804597701149)
widthsMultiGraphPave.SetY2NDC(0.85)#0.8601694915254238)

legend.SetX1NDC(0.69)
legend.SetY1NDC(0.72)
legend.SetX2NDC(0.99)
legend.SetY2NDC(0.87)
rt.gPad.Modified()

widthsMultiGraphPave.Draw()
compareRatiosCanvas.Draw()
compareRatiosCanvas.SaveAs('away-over-near.pdf')