# fits indicate not true gaussian

import ROOT as rt
import numpy as np
import matplotlib.pyplot as plt

from scipy.special import polygamma

rt.gStyle.SetOptFit(111) # show fit panel w/ fit prob, chi2/ndf, fit params, errors
rt.gStyle.SetOptStat(0)

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

######################
## Widths across fits
######################
rt.gStyle.SetLineScalePS(1)

# fit to regular gaussian
distributions = [hh08, hh12, hh20, hl08, hl12, hl20]

c_gaus = rt.TCanvas('gaus_canvas', 'gaus_canvas')
#c_gaus.SetCanvasSize(2400, 800)
c_gaus.SetCanvasSize(800, 1200)
c_gaus.Divide(2, 3)

c_gen_gaus = rt.TCanvas('gaus_gen_canvas', 'gaus_gen_canvas')
c_gen_gaus.SetCanvasSize(800, 1200)
c_gen_gaus.Divide(2, 3)

c_mises = rt.TCanvas('mises_canvas', 'gaus_mises_canvas')
c_mises.SetCanvasSize(800, 1200)
c_mises.Divide(2, 3)

c_fit_example = rt.TCanvas('fit_example', 'fit_example')
c_fit_example.SetCanvasSize(1200, 400)
c_fit_example.Divide(3, 1)

fit_functions=[FitFunction.double_gaussian, FitFunction.double_generalized_gaussian, FitFunction.double_von_mises]

gaus_params = []
gen_gaus_params = []
mises_params = []

pi = rt.TMath.Pi()
gaus_guesses =[[20000, 0.0, 1, 8000, pi, 1], [20000, 0.0, 1, 15000, pi, 1], [20000, 0.0, 1, 20000, pi, 1],
                [100, 0.0, 1, 150, pi, 1], [100, 0.0, 0.19, 250, pi, 1], [100, 0.0, 1, 400, pi, 1]]

gen_gaus_guesses = [[0.0, 0.7, 1.8, 1, pi, 1, 1.8, 1], [0.0, 0.7, 1.8, 1, pi, 1, 1.8, 1], [0.0, 0.7, 1.8, 1, pi, 1, 1.8, 1], [0.0, 0.7, 1.8, 1, pi, 1, 1.8, 1],
                        [0.0, 0.7, 1.8, 1, pi, 1, 1.8, 1], [0.0, 0.7, 1.8, 1, pi, 1, 1.8, 1], [0.0, 0.7, 1.8, 1, pi, 1, 1.8, 1], [0.0, 0.7, 1.8, 1, pi, 1, 1.8, 1]]

mises_guesses = [[20000, 0.0, 1, 8000, pi, 1], [20000, 0.0, 1, 15000, pi, 1], [20000, 0.0, 30, 20000, pi, 1],
                [100, 0.0, 1, 150, pi, 1], [100, 0.0, 0.19, 250, pi, 1], [100, 0.0, 1, 400, pi, 1]]

for i, (sparse, guess) in enumerate(zip(distributions, gaus_guesses)):
        c_gaus.cd(i + 1) 
        rt.gPad.SetLeftMargin(0.2)

        gaussian = CurveFit(sparse=sparse, dphi=sparse.make_delta_phi(), fit_function=FitFunction.double_gaussian, fit_range=FitRange.full_range, fit_params=guess)
        dphi, params = gaussian.fit()
        
        gaus_params.append(params)
        dphi.GetYaxis().SetTitle('counts')
        dphi.GetYaxis().SetTitleOffset(2.0)

        if i < 3: 
            dphi.SetMarkerColor(rt.kBlue)
            dphi.SetLineColor(rt.kBlue)
            dphi.SetTitle('Hadron-Hadron Correlation')

        else:
            dphi.SetMarkerColor(rt.kRed)
            dphi.SetLineColor(rt.kRed)
            dphi.SetTitle('Hadron-#Lambda Correlation')
        
        dphi.Draw()

        label_y_start = 0.85 #0.96 
        label_x_start = 0.3 #0.65
        pt_label_x_start = 0.07
        label_text_space = 0.06
        data_label = rt.TLatex()
        data_label.SetNDC()
        data_label.SetTextSize(0.05)
        data_label.SetTextAlign(13)
        data_label.DrawLatex(label_x_start, label_y_start, '|#eta^{trig}| < ' + str(np.abs(sparse.etaTrig[0])))
        data_label.DrawLatex(label_x_start, label_y_start - label_text_space, '|#eta^{assoc}| < ' + str(np.abs(sparse.etaAssoc[0])))
        # data_label.DrawLatex(label_x_start, label_y_start, 'pp, PYTHIA6, #sqrt{s}=7 TeV')
        # data_label.DrawLatex(label_x_start, label_y_start - label_text_space, '4 < p_{T}^{trig} < 8 GeV/c, 2 < p_{T}^{assoc} < 4 GeV/c')
        # data_label.DrawLatex(label_x_start, label_y_start - 2*label_text_space, '|#eta^{trig}| < ' + str(np.abs(sparse.etaTrig[0])) + ', |#eta^{assoc}| < ' + str(np.abs(sparse.etaAssoc[0])))

        if i == 0:
            c_fit_example.cd(1)
            rt.gPad.SetLeftMargin(0.2)
            dphi.SetTitle('Hadron-Hadron')
            dphi.GetYaxis().SetTitleOffset(2.1)
            #dphi.Draw()

            rp = rt.TRatioPlot(dphi)
            rp.SetLeftMargin(0.2)
            rp.Draw()

            label = rt.TLatex()
            label_y_start = 0.85 #0.96 
            label_x_start = 0.25
            data_label.SetTextSize(0.05)
            data_label.SetTextAlign(13)
            data_label.DrawLatex(label_x_start, label_y_start, 'Gaussian')
            data_label.DrawLatex(0.6, 0.6, '|#eta^{trig}| < ' + str(np.abs(sparse.etaTrig[0])))
            data_label.DrawLatex(0.6, 0.6 - label_text_space, '|#eta^{assoc}| < ' + str(np.abs(sparse.etaAssoc[0])))

c_gaus.Draw()

for i, (sparse, guess) in enumerate(zip(distributions, gen_gaus_guesses)):
        c_gen_gaus.cd(i + 1) 
        rt.gPad.SetLeftMargin(0.2)

        gen_gaussian = CurveFit(sparse=sparse, dphi=sparse.make_delta_phi(), fit_function=FitFunction.double_generalized_gaussian, fit_range=FitRange.full_range, fit_params=guess)
        dphi, params = gen_gaussian.fit()

        gen_gaus_params.append(params)
        dphi.GetYaxis().SetTitle('counts')
        dphi.GetYaxis().SetTitleOffset(2.0)
        dphi.Draw()

        if i < 3: 
            dphi.SetMarkerColor(rt.kBlue)
            dphi.SetLineColor(rt.kBlue)
            dphi.SetTitle('Hadron-Hadron Correlation')

        else:
            dphi.SetMarkerColor(rt.kRed)
            dphi.SetLineColor(rt.kRed)
            dphi.SetTitle('Hadron-#Lambda Correlation')
        
        dphi.Draw()

        label_y_start = 0.85 #0.96 
        label_x_start = 0.3 #0.65
        pt_label_x_start = 0.07
        label_text_space = 0.06
        data_label = rt.TLatex()
        data_label.SetNDC()
        data_label.SetTextSize(0.05)
        data_label.SetTextAlign(13)
        data_label.DrawLatex(label_x_start, label_y_start, '|#eta^{trig}| < ' + str(np.abs(sparse.etaTrig[0])))
        data_label.DrawLatex(label_x_start, label_y_start - label_text_space, '|#eta^{assoc}| < ' + str(np.abs(sparse.etaAssoc[0])))
        # data_label.DrawLatex(label_x_start, label_y_start, 'pp, PYTHIA6, #sqrt{s}=7 TeV')
        # data_label.DrawLatex(label_x_start, label_y_start - label_text_space, '4 < p_{T}^{trig} < 8 GeV/c, 2 < p_{T}^{assoc} < 4 GeV/c')
        # data_label.DrawLatex(label_x_start, label_y_start - 2*label_text_space, '|#eta^{trig}| < ' + str(np.abs(sparse.etaTrig[0])) + ', |#eta^{assoc}| < ' + str(np.abs(sparse.etaAssoc[0])))

        if i == 0:
            c_fit_example.cd(2)
            rt.gPad.SetLeftMargin(0.2)
            dphi.GetYaxis().SetTitleOffset(2.1)

            #dphi.Draw()
            rp2 = rt.TRatioPlot(dphi)
            rp2.SetLeftMargin(0.2)
            rp2.Draw()

            label = rt.TLatex()
            label_y_start = 0.85 #0.96 
            label_x_start = 0.25
            data_label.SetTextSize(0.05)
            data_label.SetTextAlign(13)
            data_label.DrawLatex(label_x_start, label_y_start, 'Gen. Gaussian')
            data_label.DrawLatex(0.6, 0.6, '|#eta^{trig}| < ' + str(np.abs(sparse.etaTrig[0])))
            data_label.DrawLatex(0.6, 0.6 - label_text_space, '|#eta^{assoc}| < ' + str(np.abs(sparse.etaAssoc[0])))


for i, (sparse, guess) in enumerate(zip(distributions, mises_guesses)):
        c_mises.cd(i + 1) 
        rt.gPad.SetLeftMargin(0.2)

        mises = CurveFit(sparse=sparse, dphi=sparse.make_delta_phi(), fit_function=FitFunction.double_von_mises, fit_range=FitRange.full_range, fit_params=guess)
        dphi, params = mises.fit()
        mises_params.append(params)

        dphi.GetYaxis().SetTitle('counts')
        dphi.GetYaxis().SetTitleOffset(2.0)
        dphi.Draw()

        if i < 3: 
            dphi.SetMarkerColor(rt.kBlue)
            dphi.SetLineColor(rt.kBlue)
            dphi.SetTitle('Hadron-Hadron Correlation')

        else:
            dphi.SetMarkerColor(rt.kRed)
            dphi.SetLineColor(rt.kRed)
            dphi.SetTitle('Hadron-#Lambda Correlation')
        
        dphi.Draw()

        label_y_start = 0.85 #0.96 
        label_x_start = 0.3 #0.65
        pt_label_x_start = 0.07
        label_text_space = 0.06
        data_label = rt.TLatex()
        data_label.SetNDC()
        data_label.SetTextSize(0.05)
        data_label.SetTextAlign(13)
        data_label.DrawLatex(label_x_start, label_y_start, '|#eta^{trig}| < ' + str(np.abs(sparse.etaTrig[0])))
        data_label.DrawLatex(label_x_start, label_y_start - label_text_space, '|#eta^{assoc}| < ' + str(np.abs(sparse.etaAssoc[0])))
        # data_label.DrawLatex(label_x_start, label_y_start, 'pp, PYTHIA6, #sqrt{s}=7 TeV')
        # data_label.DrawLatex(label_x_start, label_y_start - label_text_space, '4 < p_{T}^{trig} < 8 GeV/c, 2 < p_{T}^{assoc} < 4 GeV/c')
        # data_label.DrawLatex(label_x_start, label_y_start - 2*label_text_space, '|#eta^{trig}| < ' + str(np.abs(sparse.etaTrig[0])) + ', |#eta^{assoc}| < ' + str(np.abs(sparse.etaAssoc[0])))

        if i == 0:
            c_fit_example.cd(3)
            rt.gPad.SetLeftMargin(0.2)
            dphi.GetYaxis().SetTitleOffset(2.1)

            #dphi.Draw()
            rp3 = rt.TRatioPlot(dphi)
            rp3.SetLeftMargin(0.2)
            rp3.Draw()

            label = rt.TLatex()
            label_y_start = 0.85 #0.96 
            label_x_start = 0.25
            data_label.SetTextSize(0.05)
            data_label.SetTextAlign(13)
            label = rt.TLatex()
            data_label.DrawLatex(label_x_start, label_y_start, 'Von Mises')
            data_label.DrawLatex(0.6, 0.6, '|#eta^{trig}| < ' + str(np.abs(sparse.etaTrig[0])))
            data_label.DrawLatex(0.6, 0.6 - label_text_space, '|#eta^{assoc}| < ' + str(np.abs(sparse.etaAssoc[0])))
        

c_gaus.SaveAs('gaussian_fit.pdf')
c_gen_gaus.SaveAs('gen_gaussian_fit.pdf')
c_mises.SaveAs('mises_fit.pdf')

######################
## Width Calculations
######################

gaus_params = np.array(gaus_params)
gaus_stdev_near = gaus_params[:, 2, 0]
gaus_stdev_away = gaus_params[:, 5, 0]

gaus_stdev_near_err = gaus_params[:, 2, 1]
gaus_stdev_away_err = gaus_params[:, 5, 1]


def gen_gaus_stdev(alpha, beta):
    '''
    Returns the standard deviation of a generalized gaussian, given alpha and beta. 
    '''
    var = alpha**2 * rt.Math.tgamma(3 / beta) / rt.Math.tgamma(1 / beta)
    
    return rt.TMath.Sqrt(var)

def gen_gaus_stdev_err(alpha, alpha_err, beta, beta_err):
     '''
     Returns the error in the standard deviation of a generalized gaussain, given alpha and beta---and their uncertainties.

     Partial derivatives calculated with mathematica. 
     '''

     stdev = gen_gaus_stdev(alpha, beta)

     partial_alpha = stdev / alpha

     partial_beta_num = stdev * (polygamma(0, 1/beta) - 3 * polygamma(0, 3/beta))
     partial_beta_denom = 2 * beta**2
     partial_beta = partial_beta_num / partial_beta_denom

     err_squared = (partial_alpha * alpha_err)**2 + (partial_beta * beta_err)**2 

     return np.sqrt(err_squared)


gen_gaus_params = np.array(gen_gaus_params)

gen_gaus_params_alpha_near = gen_gaus_params[:, 1, 0]
gen_gaus_params_beta_near = gen_gaus_params[:, 2, 0] 

gen_gaus_params_alpha_away = gen_gaus_params[:, 5, 0]
gen_gaus_params_beta_away = gen_gaus_params[:, 6, 0]

gen_gaus_params_alpha_near_err = gen_gaus_params[:, 1, 1]
gen_gaus_params_beta_near_err = gen_gaus_params[:, 2, 1] 

gen_gaus_params_alpha_away_err = gen_gaus_params[:, 5, 1]
gen_gaus_params_beta_away_err = gen_gaus_params[:, 6, 1]

gen_gaus_stdev_near = np.array([gen_gaus_stdev(alpha, beta) for (alpha, beta) in zip(gen_gaus_params_alpha_near, gen_gaus_params_beta_near)])
gen_gaus_stdev_away = np.array([gen_gaus_stdev(alpha, beta) for (alpha, beta) in zip(gen_gaus_params_alpha_away, gen_gaus_params_beta_away)])

gen_gaus_stdev_near_err = np.array([gen_gaus_stdev_err(alpha, alpha_err, beta, beta_err) for (alpha, alpha_err, beta, beta_err) in zip(gen_gaus_params_alpha_near, gen_gaus_params_alpha_near_err, gen_gaus_params_beta_near, gen_gaus_params_beta_near_err)])
gen_gaus_stdev_away_err = np.array([gen_gaus_stdev_err(alpha, alpha_err, beta, beta_err) for (alpha, alpha_err, beta, beta_err) in zip(gen_gaus_params_alpha_away, gen_gaus_params_alpha_away_err, gen_gaus_params_beta_away, gen_gaus_params_beta_away_err)])

def mises_stdev(kappa):
    '''
    Returns the standard deviation of a von mises distribution, given alpha and beta.
    '''
    var = 1 - rt.TMath.BesselI1(kappa) / rt.TMath.BesselI0(kappa)
    var *= 2
    #var = -2 * rt.TMath.Log(rt.TMath.BesselI1(kappa) / rt.TMath.BesselI0(kappa))

    return rt.TMath.Sqrt(var)

def mises_stdev_err(kappa, kappa_err):
     '''
     Returns the error in the standard deviation of the Von Mises distribution, given kappa and its error. 
     
     This is used by ALICE in published results. Since the published formula returns negative values, we take the absolute value. 
     '''
     I1 = rt.TMath.BesselI1(kappa)
     I0 = rt.TMath.BesselI0(kappa)
     
     stdev = mises_stdev(kappa)
     err = (1 / stdev) * (I1/I0 - I0/I1 + 1/kappa) * kappa_err

     return np.abs(err)

# I1NS = TMath::BesselI1(NS_Kappa);
# I0NS = TMath::BesselI0(NS_Kappa);
# NSSigma = TMath::Sqrt(-2*TMath::Log(I1NS/I0NS));
# NSSigError = (1./NSSigma)*(I1NS/I0NS - I0NS/I1NS + 1./NS_Kappa)*(NS_Kappa);

mises_params = np.array(mises_params)
mises_kappa_near = mises_params[:, 2, 0]
mises_kappa_away = mises_params[:, 5, 0]

mises_kappa_near_err = mises_params[:, 2, 1]
mises_kappa_away_err = mises_params[:, 5, 1]

mises_stdev_near = np.array([mises_stdev(kappa) for kappa in mises_kappa_near])
mises_stdev_away = np.array([mises_stdev(kappa) for kappa in mises_kappa_away])

mises_stdev_near_err = np.array([mises_stdev_err(kappa, kappa_err) for (kappa, kappa_err) in zip(mises_kappa_near, mises_kappa_near_err)])
mises_stdev_away_err = np.array([mises_stdev_err(kappa, kappa_err) for (kappa, kappa_err) in zip(mises_kappa_away, mises_kappa_away_err)])

print('near', mises_stdev_near)
print('near err', mises_stdev_near_err)

SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 12

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

plt.rcParams['figure.dpi'] = 400

fig = plt.figure()
gs = fig.add_gridspec(2, 2, hspace=0.2, wspace=0.3)
((ax1, ax2), (ax3, ax4)) = gs.subplots(sharex='col')#, sharey='row')

eta = np.array([0.8, 1.2, 2.0])
offset = 0.08
#ax1.scatter(eta, gaus_stdev_near[:3], color='k')
#ax1.scatter(eta, gen_gaus_stdev_near[:3], color='b')
#ax1.scatter(eta, mises_stdev_near[:3], color='r')

ax1.errorbar(eta - offset, gaus_stdev_near[:3], yerr=gaus_stdev_near_err[:3], fmt='sk')#, color='k')
ax1.errorbar(eta, gen_gaus_stdev_near[:3], yerr=gen_gaus_stdev_near_err[:3], fmt='sb')#, color='b')
ax1.errorbar(eta + offset, mises_stdev_near[:3], mises_stdev_near_err[:3], fmt='sr')#, color='r')

#ax2.scatter(eta, gaus_stdev_near[3:], color='k', label='Gaussian')
#ax2.scatter(eta, gen_gaus_stdev_near[3:], color='b', label='Gen. Gaussian')
#ax2.scatter(eta, mises_stdev_near[3:], color='r',  label='Von Mises')

ax2.errorbar(eta - offset, gaus_stdev_near[3:], yerr=gaus_stdev_near_err[3:], fmt='sk', label='Gaussian')
ax2.errorbar(eta, gen_gaus_stdev_near[3:], yerr=gen_gaus_stdev_near_err[3:], fmt='sb', label='Gen. Gaussian')
ax2.errorbar(eta + offset, mises_stdev_near[3:], yerr=mises_stdev_near_err[3:], fmt='sr', label='Von Mises')

# ax3.scatter(eta, gaus_stdev_away[:3], color='k')
# ax3.scatter(eta, gen_gaus_stdev_away[:3], color='b')
# ax3.scatter(eta, mises_stdev_away[:3], color='r')

ax3.errorbar(eta - offset, gaus_stdev_away[:3], yerr=gaus_stdev_away_err[:3], fmt='ok')
ax3.errorbar(eta, gen_gaus_stdev_away[:3], yerr=gen_gaus_stdev_away_err[:3], fmt='ob')
ax3.errorbar(eta + offset, mises_stdev_away[:3], yerr=mises_stdev_away_err[:3], fmt='or')

# ax4.scatter(eta, gaus_stdev_away[3:], color='k')
# ax4.scatter(eta, gen_gaus_stdev_away[3:], color='b')
# ax4.scatter(eta, mises_stdev_away[3:], color='r')

ax4.errorbar(eta - offset, gaus_stdev_away[3:], yerr=gaus_stdev_away_err[3:], fmt='ok')#, color='k')
ax4.errorbar(eta, gen_gaus_stdev_away[3:], yerr=gen_gaus_stdev_away_err[3:], fmt='ob')#, color='b')
ax4.errorbar(eta + offset, mises_stdev_away[3:], mises_stdev_away_err[3:], fmt='or')#, color='r')

ax3.set_xlim(0.6, 2.2)
ax4.set_xlim(0.6, 2.2)
ax3.set_xticks(eta)
ax4.set_xticks(eta)

ax1.set_ylim(0.17, 0.20)
ax2.set_ylim(0.15, 0.28)
ax3.set_ylim(0.34, 0.4)
ax4.set_ylim(0.30, 0.46)

ax1.set_ylabel('$\sigma_{NS}$')
ax3.set_ylabel('$\sigma_{AS}$')

ax3.set_xlabel('$|\eta|<x$')
ax4.set_xlabel('$|\eta|<x$')

ax1.set_title('h$-$h')
ax2.set_title('h$-\Lambda$')

ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5), markerscale=0)

#plt.tight_layout()
fig.suptitle('$\sigma_{\mathrm{NS}}$ and $\sigma_{\mathrm{AS}}$ Across Cuts and Fits')
plt.savefig('widths-across-fits.pdf', bbox_inches='tight')


fig, axs = plt.subplots(2, sharex=True)

def ratio_err(numerator, numerator_err, denomenator, denomenator_err):
     ratio = numerator / denomenator
     
     return ratio * np.sqrt((numerator_err / numerator)**2 + (denomenator_err / denomenator)**2)

eta = np.array([0.8, 1.2, 2.0])
offsets = np.array([-0.05, 0.0, 0.05])
ns_stdevs = [gaus_stdev_near, gen_gaus_stdev_near, mises_stdev_near]
ns_stdevs_err = [gaus_stdev_near_err, gen_gaus_stdev_near_err, mises_stdev_near_err]
as_stdevs = [gaus_stdev_away, gen_gaus_stdev_away, mises_stdev_away]
as_stdevs_err = [gaus_stdev_away_err, gen_gaus_stdev_away_err, mises_stdev_away_err]
format = [['sk', 'ok'], ['sb', 'ob'], ['sr', 'or']]
labels = ['Gaussian', 'Gen. Gaussian', 'Von Mises']

for ns_stdev, ns_stdev_err, as_stdev, as_stdev_err, fmt, label, offset in zip(ns_stdevs, ns_stdevs_err, as_stdevs, as_stdevs_err, format, labels, offsets):
    near_ratios = ns_stdev[3:] / ns_stdev[:3] # h-L / h-h near side
    away_ratios = as_stdev[3:] / as_stdev[:3] # h-L / h-h away side

    near_ratios_err = ratio_err(ns_stdev[3:], ns_stdev_err[3:], ns_stdev[:3], ns_stdev_err[:3])#ns_stdev_err[3:] / ns_stdev_err[:3]
    away_ratios_err = ratio_err(as_stdev[3:], as_stdev_err[3:], as_stdev[:3], as_stdev_err[:3])#as_stdev_err[3:] / as_stdev_err[:3]
    print('err', near_ratios_err)
    print(away_ratios_err)
    near_fmt, away_fmt = fmt 
    axs[0].errorbar(eta + offset, near_ratios, yerr=near_ratios_err, fmt=near_fmt, label=label)
    axs[1].errorbar(eta + offset, away_ratios, yerr=away_ratios_err, fmt=away_fmt, label=label)

axs[0].set_xticks(eta)
axs[1].set_xticks(eta)

#plt.legend(markerscale=0)
axs[0].legend(loc='center left', bbox_to_anchor=(1, 0.5), markerscale=0)

axs[0].set_ylim(0.8, 1.6)
axs[1].set_ylim(0.8, 1.2)

axs[1].set_xlabel('$|\eta|<x$')
axs[0].set_ylabel('$\sigma^{h-\Lambda}_{NS} / \sigma^{h-h}_{NS}$')
axs[1].set_ylabel('$\sigma^{h-\Lambda}_{AS} / \sigma^{h-h}_{AS}$')

fig.suptitle('Width Ratios for Various Fits and $\eta$ Cuts')

#plt.tight_layout()
plt.savefig('ratios-across-fits.pdf', bbox_inches='tight')

##################
## Width Examples
##################

c_fit_example.SaveAs('example_fits.pdf')
