import ROOT as rt
import numpy as np

from enum import Enum
from formatting import process_tf1


class FitFunction(Enum):
    gaussian = 'gaussian'
    generalized_gassian = 'generalized_gaussian' 
    von_mises = 'von_mises'
    double_gaussian = 'double_gaussian'
    double_generalized_gaussian = 'double_gen_gaussian'
    double_von_mises = 'double_von_mises'


class FitRange(Enum):
    away_side = 'away_side'
    near_side = 'near_side'
    full_range = 'full_range'


class CurveFit:
    def __init__(self, sparse, dphi, fit_function, fit_range, fit_params, ratio_plot=False):
        '''
        Defines an easy-to-use API for fitting delta phi distributions for various fit functions and ranges. The internals are ugly, but it is easily maintainable.

        Attributes:
            sparse (SparseAnalyzer): takes name from SparseAnalyzer, might make this generate the dphi distribution automatically
            dphi (TH1): delta phi distribution to fit
            fit_function (Enum): option from FitFunction. Defines the fit function that gets used.
            fit_range (Enum): option from FitRange. Defines the range we want to fit over (near, away, entire domain).
            fit_params (array-like): guess parameters for fit. Guess for means are used to fix the mean in fits. 
            ratio_plot (boolean): True if want a ratio plot from the fit method. False if not. 
        
        To do:
        Think abt if should build dphi in constructor.
        '''
        self.name = sparse.name
        self.dphi = dphi
        self.f = None
        self.fit_range = fit_range
        self.ratio_plot = ratio_plot
        self.fit_name = f'{self.name}_{fit_function.value}_{fit_range.value}'

        self.dphi.SetName(f'{self.fit_name}_dphi')
        self.dphi.GetYaxis().SetRangeUser(0, 1.3 * self.dphi.GetMaximum())

        pi = rt.TMath.Pi()

        # set what range we fit over
        if fit_range is FitRange.full_range:
            self.fit_range = [-pi/2, 3*pi/2]

        elif fit_range is FitRange.near_side:
            self.fit_range = [-pi/2, pi/2]

        elif fit_range is FitRange.away_side:
            self.fit_range = [pi/2, 3*pi/2]
        
        # define fit function, fix params and set param ranges
        if fit_function is FitFunction.gaussian:

            self.f = rt.TF1(f'{self.fit_name}_fit', 'gaus', *self.fit_range)
            self._num_params = 3

            amp, mean, stdev = fit_params
            self.f.SetParameters(amp, mean, stdev)
            self.f.FixParameter(1, mean) # fix mean
            self.f.SetParLimits(2, 0.0, 100) # ensure stdev is positive to reduce param space

        elif fit_function is FitFunction.generalized_gassian:

            # 0: mean, 1: alpha, 2: beta, 3: amp
            self.f = rt.TF1(f'{self.fit_name}_fit', "[3]*([2]/(2*[1]*TMath::Gamma(1/[2])))*TMath::Exp(-TMath::Power(TMath::Abs(x-[0])/[1],[2]))",
                            *self.fit_range) 
            self._num_params = 4           

            mean, alpha, beta, amp = fit_params
            self.f.SetParameters(mean, alpha, beta, amp)
            self.f.FixParameter(0, mean) # fix mean 
            self.f.SetParLimits(2, 1, 2) # limit beta
        
        elif fit_function is FitFunction.von_mises:
            
            # 0: Amp/Yield, 1: mean, 2: kappa
            self.f = rt.TF1(f"{self.fit_name}_fit", "[0]/(2*TMath::Pi()*TMath::BesselI0([2]))*TMath::Exp([2]*TMath::Cos(x- 2*TMath::Pi() - [1]))",
                        *self.fit_range)
            self._num_params = 3
            
            amp, mean, kappa = fit_params
            self.f.SetParameters(amp, mean, kappa)
            self.f.FixParameter(1, mean) # fix mean

        elif fit_function is FitFunction.double_gaussian:
            
            self.f = rt.TF1(f'{self.fit_name}_fit', 'gaus(0) + gaus(3)', *self.fit_range)
            self._num_params = 6

            amp_near, mean_near, stdev_near, amp_away, mean_away, stdev_away = fit_params
            self.f.SetParameters(amp_near, mean_near, stdev_near, 
                                    amp_away, mean_away, stdev_away)
            self.f.FixParameter(1, mean_near)
            self.f.FixParameter(4, mean_away)
            self.f.SetParLimits(2, 0, 100)
            self.f.SetParLimits(5, 0, 100)

        elif fit_function is FitFunction.double_generalized_gaussian:
            
            # 0: Mean NS, 1: alpha NS, 2: beta NS, 3: Amp GGaus NS, 4: Mean AS, 5: alpha AS, 6: beta AS, 7: Amp GGaus AS 
            self.f = rt.TF1(f'{self.fit_name}_fit', "[3]*([2]/(2*[1]*TMath::Gamma(1/[2])))*TMath::Exp(-TMath::Power(TMath::Abs(x-[0])/[1],[2])) + [7]*([6]/(2*[5]*TMath::Gamma(1/[6])))*TMath::Exp(-TMath::Power(TMath::Abs(x-[4])/[5],[6]))",
                        *self.fit_range)
            self._num_params = 8

            mean_near, alpha_near, beta_near, amp_near, mean_away, alpha_away, beta_away, amp_away = fit_params
            self.f.SetParameters(mean_near, alpha_near, beta_near, amp_near, mean_away, alpha_away, beta_away, amp_away) 

            self.f.FixParameter(0, mean_near)
            self.f.FixParameter(4, mean_away)

            self.f.SetParLimits(2, 1, 2)
            self.f.SetParLimits(6, 1, 2)

        elif fit_function is FitFunction.double_von_mises:
            
            # 0: NS Amp/Yield, 1: NS mean, 2: NS kappa, 3: AS Amp/Yield, 4: AS mean, 5: AS kappa
            self.f = rt.TF1(f'{self.fit_name}_fit', "[0]/(2*TMath::Pi()*TMath::BesselI0([2]))*TMath::Exp([2]*TMath::Cos(x- 2*TMath::Pi() - [1])) + [3]/(2*TMath::Pi()*TMath::BesselI0([5]))*TMath::Exp([5]*TMath::Cos(x- 2*TMath::Pi()-[4]))",
                        *self.fit_range)
            self._num_params = 6

            amp_near, mean_near, kappa_near, amp_away, mean_away, kappa_away = fit_params
            self.f.SetParameters(amp_near, mean_near, kappa_near, amp_away, mean_away, kappa_away)
            # Fix parameter 1 at 0; fix 4 at pi
            self.f.FixParameter(1, mean_near)
            self.f.FixParameter(4, mean_away)            
        
        process_tf1(self.f)
        self.f.SetNpx(1000)

        
    def fit(self):
        '''
        Apply fit to delta phi distribution.

        Returns
            Fitted TH1 or TRatioPlot: TRatioPlot if ratio_plot is true, TH1 if not.
            params_w_err (numpy array): Parameters and errors, shape=(number of parameters, 2)
        '''
        self.dphi.Fit(f'{self.fit_name}_fit', 'BR')

        params = [self.f.GetParameter(i) for i in range(self._num_params)]
        errs = [self.f.GetParError(i) for i in range(self._num_params)]
        params_w_err = np.array(list(zip(params, errs)), dtype='d')

        if self.ratio_plot:
            return rt.TRatioPlot(self.dphi), params_w_err
        else:
            return self.dphi, params_w_err