from enum import Enum

eps = 0.0001
class SparseAnalyzer:

    def __init__(self, sparse, name, etaAssoc=[-0.8, 0.8-eps], etaTrig=[-0.8, 0.8-eps], ptAssoc=[2, 4-eps], ptTrig=[4, 8-eps], apply_cut=True):
        '''
        Bundles together THnSparse with automatically applied pT and eta cuts, plus useful methods for analysis. Important: apply_cut only toggles the trigger and associate eta cuts. 
        Each method first calls a cutting function, and finishes off with a function that zeros out these cuts. This ensures we can treat our THnSparse like a multidimensional 
        histogram with all cuts already applied, without worrying about changing the underlying THnSparse that ROOT stores a reference to. 

        Attributes
            sparse (THnSparse): THnSparse with 6 axes ordered as: trigger eta, associate eta, trigger pT, associate pT, delta phi, delta eta
            name (string): Name used to generate ROOT names for our THnSparse and histograms. This name should be descriptive. 
                For example, for hadron-hadron data with an eta=[-0.8, 0.8] cut, we might use the name 'hh08'.
            etaAssoc (array-like): Associate single-particle eta cut. Subtract a small epsilon from the right edge to avoid landing exactly on a right bin edge.
            etaTrig (array-like): Trigger single-particle eta cut. Subtract a small epsilon from the right edge to avoid landing exactly on a right bin edge.
            ptAssoc (array-like): Associate single-particle pT cut. Subtract a small epsilon from the right edge to avoid landing exactly on a right bin edge.
            ptTrig (array-like): Trigger single-particle pT cut. Subtract a small epsilon from the right edge to avoid landing exactly on a right bin edge.
            apply_cut (bool): True if we desire a cut on eta, False if not. Note: this does not affect the pT cut.
        '''
        self.sparse = sparse
        self.name = name
        self.apply_cut = apply_cut
        self._axis = {'trigger eta': 0, 'associate eta': 1, 'trigger pT': 2, 'associate pT': 3, 'delta phi': 4, 'delta eta': 5}


        self.ptTrig = ptTrig
        self.ptAssoc = ptAssoc

        if self.apply_cut:
            self.etaTrig = etaTrig
            self.etaAssoc = etaAssoc
            self._cut()

        
    def _cut(self):
        '''
        Applies cut to self.sparse.
        ''' 
        ptTrigMin, ptTrigMax = self.ptTrig
        ptAssocMin, ptAssocMax = self.ptAssoc

        self.sparse.GetAxis(self._axis['trigger pT']).SetRangeUser(ptTrigMin, ptTrigMax)
        self.sparse.GetAxis(self._axis['associate pT']).SetRangeUser(ptAssocMin, ptAssocMax)
        
        if not self.apply_cut:
            return None
        
        etaTrigMin, etaTrigMax = self.etaTrig
        etaAssocMin, etaAssocMax = self.etaAssoc 

        self.sparse.GetAxis(self._axis['trigger eta']).SetRangeUser(etaTrigMin, etaTrigMax)
        self.sparse.GetAxis(self._axis['associate eta']).SetRangeUser(etaAssocMin, etaAssocMax)

        return None


    def _reset_axes(self):    
        '''
        Resets cut on self.sparse.
        ''' 
        self.sparse.GetAxis(self._axis['trigger pT']).SetRangeUser(0, 0)
        self.sparse.GetAxis(self._axis['associate pT']).SetRangeUser(0, 0)
        self.sparse.GetAxis(self._axis['trigger eta']).SetRangeUser(0, 0)
        self.sparse.GetAxis(self._axis['associate eta']).SetRangeUser(0, 0)


    @staticmethod    
    def _format_dist(dist, color=1, size=1.4, style=20):
        '''
        Format histogram for readability. 
        For configuring color, size, style, see: https://root.cern.ch/doc/master/classTAttMarker.html
        
        Args
            dist (TH1, TH2): Histogram to format.
            color (int): Set marker color.
            size (float): Set marker size
            style (int): Set marker style. 
        '''

        #gPad.SetTickx();
        #gPad.SetTicky();
        dist.SetMarkerSize(size)
        dist.SetMarkerColor(color)
        dist.SetLineColor(color)
        dist.SetMarkerStyle(style)

        dist.GetXaxis().SetTitleFont(42)
        dist.GetXaxis().SetTitleOffset(1.0)
        dist.GetXaxis().SetTitleSize(0.045)
        dist.GetXaxis().SetLabelSize(0.045)

        dist.GetYaxis().SetTitleOffset(1.1)
        dist.GetYaxis().SetTitleSize(0.05)
        dist.GetYaxis().SetLabelSize(0.05)
        dist.GetYaxis().SetLabelFont(42)
        dist.GetXaxis().SetLabelFont(42)
        dist.GetYaxis().SetTitleFont(42)
        

    @staticmethod
    def _process_pave(pave, font=42, size=0.06):
        pave.SetFillStyle(0)
        pave.SetBorderSize(0)
        #pave.SetTextColor(kBlack)
        pave.SetTextSize(size)
        pave.SetTextFont(font)

        return pave


    def make_pave(self):
        pave = rt.TPaveText(0.1, 0.75, 0.5, 0.9, 'NDC')
        pave = self._process_pave(pave)

        pave.AddText('pp, PYTHIA6, #sqrt{s}=14 TeV')
        pave.AddText('4 < p_{T}^{trig} < 8 GeV/c, 2 < p_{T}^{assoc} < 4 GeV/c')
    
        if self.apply_cut:
            pave.AddText('|#eta^{trig}| < ' + str(np.abs(self.etaTrig[0])) + ', |#eta^{assoc}| < ' + str(np.abs(self.etaAssoc[0])))
    
        return pave


    def make_1D_hist(self, axis):
        '''
        Make a 1D distribution from a given axis of self.sparse.

        Args
            axis (string): Name of axis to project. Must be 

        Returns
        '''
        self._cut()
        name = {'trigger eta': '_trig_eta', 'associate eta': '_assoc_eta', 'trigger pT': '_trig_pt', 'associate pT': '_assoc_pt', 
                'delta phi': '_dphi', 'delta eta': '_deta'}

        hist = self.sparse.Projection(self._axis[axis])
        hist.GetYaxis().SetRangeUser(0, 1.3 * hist.GetMaximum())
        hist.SetName(f'{self.name}{name[axis]}')
        self._format_dist(hist)
        
        self._reset_axes()

        return hist


    def make_delta_eta(self):
        self._cut()
        dEta = self.sparse.Projection(self._axis['delta eta'])
        dEta.GetYaxis().SetRangeUser(0, 1.3 * dEta.GetMaximum())
        dEta.SetName(f'{self.name}DEta')
        self._format_dist(dEta)

        self._reset_axes()

        return dEta
    
    
    def make_delta_phi(self):
        self._cut()        

        dPhi = self.sparse.Projection(self._axis['delta phi'])
        
        dPhi.GetYaxis().SetRangeUser(0, 1.3 * dPhi.GetMaximum())
        self._format_dist(dPhi)
        dPhi.SetName(f'{self.name}DPhi')
        dPhi.GetXaxis().SetTitle('#Delta#varphi (rad)')

        self._reset_axes()

        return dPhi
    
        # raise NotImplementedError, 'implement pave'
    
    
    def make_2D_correlation(self):
        self._cut()        
        dPhiDEta = self.sparse.Projection(self._axis['delta eta'], self._axis['delta phi'])
        dPhiDEta.GetYaxis().SetRangeUser(-1.2, 1.2-eps)
        self._format_dist(dPhiDEta)
        dPhiDEta.SetName(f'{self.name}DPhiDEta')

        self._reset_axes()

        return dPhiDEta
    
    
    def integrate_away_side(self):
        '''
        Helper function for get_away_side_ratios(). Calculates the counts in the away side peak of a 2D delta phi delta eta correlation.

        Args
            correlation (TH2):

        Returns
            away_yield (float): 
        '''
        # bounds of integral
        pi = rt.TMath.Pi()
        eps = 0.0001 # in case upper bound lands on right bin edge
        dPhi_lower = pi / 2
        dPhi_upper = 3 * pi / 2 - eps
        dEta_lower = -50
        dEta_upper = 50 - eps

        correlation = self.make_2D_correlation()

        dPhiLowerBin = correlation.GetXaxis().FindBin(dPhi_lower)
        dPhiUpperBin = correlation.GetXaxis().FindBin(dPhi_upper)

        dEtaLowerBin = correlation.GetYaxis().FindBin(dEta_lower)
        dEtaUpperBin = correlation.GetYaxis().FindBin(dEta_upper)

        away_yield = correlation.Integral(dPhiLowerBin, dPhiUpperBin, dEtaLowerBin, dEtaUpperBin)

        return away_yield
    
    
    def integrate_near_side(self):
        '''
        Calculates the counts in the near side peak of a 2D delta phi delta eta correlation.

        Args
            correlation (TH2):

        Returns
            near_yield (float): 
        '''
        # bounds of integral
        pi = rt.TMath.Pi()
        eps = 0.0001 # in case upper bound lands on right bin edge
        dPhi_lower = -pi / 2
        dPhi_upper = pi / 2 - eps
        dEta_lower = -50
        dEta_upper = 50 - eps

        correlation = self.make_2D_correlation()

        dPhiLowerBin = correlation.GetXaxis().FindBin(dPhi_lower)
        dPhiUpperBin = correlation.GetXaxis().FindBin(dPhi_upper)

        dEtaLowerBin = correlation.GetYaxis().FindBin(dEta_lower)
        dEtaUpperBin = correlation.GetYaxis().FindBin(dEta_upper)

        near_yield = correlation.Integral(dPhiLowerBin, dPhiUpperBin, dEtaLowerBin, dEtaUpperBin)

        return near_yield
    
    
    def fit_delta_phi(self, dPhi, fit_params, ratio_plot=False):
        '''
        Args
            dPhi (TH1): delta phi distribution
            fit_params (array_like): Guess parameters. The order of the parameters is: first gaussian amplitude, first gaussian stdev,
                second gaussian amplitude, second gaussian stdev.
            ratio_plot (bool): True if desire a ratio plot, False if only want to fit the TH1. Use this option if you want to analyze residuals. 

        Returns
            None or ratio plot: if ratio_plot is True, returns a ROOT ratio plot. Else, returns None. 
        '''

        dPhi.GetYaxis().SetRangeUser(0, 1.3 * dPhi.GetMaximum())
        #self._format_dist(dPhi)
        process_hist(dPhi)
        dPhi.SetName(f'{self.name}_dPhi_fitted')
        
        doubleGauss = rt.TF1(f'{self.name}_fit', 'gaus(0) + gaus(3)')
        process_tf1(doubleGauss)
        doubleGauss.SetNpx(1000)

        amp_near, stdev_near, amp_away, stdev_away = fit_params
        doubleGauss.SetParameters(amp_near, 0.0, stdev_near, 
                                  amp_away, rt.TMath.Pi(), stdev_away)
        doubleGauss.FixParameter(1, 0.0)
        doubleGauss.FixParameter(4, rt.TMath.Pi())
        # doubleGauss.SetParLimits(2, 0, 100)
        # doubleGauss.SetParLimits(5, 0, 100)
        #doubleGauss.FixParameter(0.18)
        doubleGauss.SetParLimits(5, 0, 100)

        dPhi.Fit(f'{self.name}_fit', 'B')

        
        fit = dPhi.GetFunction(f'{self.name}_fit')
        params = np.array(list(zip([fit.GetParameter(i) for i in range(6)], [fit.GetParError(i) for i in range(6)])), dtype='d')
        
        if ratio_plot:
            # rp.GetLowerRefYaxis().SetTitle("ratio")
            #rp.GetUpperRefYaxis().SetTitle("entries") # kills kernel for unknown reason, maybe because projected from sparse?
            # raise NotImplementedError('need to make sure rp works w/ params')
            return rt.TRatioPlot(dPhi), params
                
        else:
            return dPhi, params
    
    
    def fit_von_mises(self, dPhi, fit_params, ratio_plot=False):
        '''
        Same num of params as gaussian.
        '''
        dPhi.SetName(f'{self.name}_dPhi_fitted_mises')

        # 0: NS Amp/Yield, 1: NS mean, 2: NS kappa, 3: AS Amp/Yield, 4: AS mean, 5: AS kappa
        mises = rt.TF1(f"{self.name}_von_mises", "[0]/(2*TMath::Pi()*TMath::BesselI0([2]))*TMath::Exp([2]*TMath::Cos(x- 2*TMath::Pi() - [1])) + [3]/(2*TMath::Pi()*TMath::BesselI0([5]))*TMath::Exp([5]*TMath::Cos(x- 2*TMath::Pi()-[4]))",
                       -rt.TMath.Pi()/2,
                       +(2)*rt.TMath.Pi())
        process_tf1(mises)
        mises.SetNpx(1000)

        amp_near, kappa_near, amp_away, kappa_away = fit_params
        mises.SetParameters(amp_near, 0.0, kappa_near, amp_away, rt.TMath.Pi(), kappa_away)
        # Fix parameter 1 at 0; fix 4 at pi
        mises.FixParameter(1, 0.0)
        mises.FixParameter(4, rt.TMath.Pi())

        dPhi.Fit(f'{self.name}_von_mises', 'B')

        fit = dPhi.GetFunction(f'{self.name}_von_mises')
        params = np.array(list(zip([fit.GetParameter(i) for i in range(6)], [fit.GetParError(i) for i in range(6)])), dtype='d')  
    
        if ratio_plot:
            return rt.TRatioPlot(dPhi), params
                
        else:
            return dPhi, params
    

    def fit_near_von_mises(self, dPhi, fit_params, ratio_plot=False):
        '''
        Fit near-side with Von Mises. 
        '''
        dPhi.SetName(f'{self.name}_dPhi_fitted_near_mises')

        # 0: NS Amp/Yield, 1: NS mean, 2: NS kappa
        mises = rt.TF1(f"{self.name}_von_mises", "[0]/(2*TMath::Pi()*TMath::BesselI0([2]))*TMath::Exp([2]*TMath::Cos(x- 2*TMath::Pi() - [1]))",
                       -rt.TMath.Pi()/2,
                       rt.TMath.Pi()/2)
        process_tf1(mises)
        mises.SetNpx(1000)

        amp_near, kappa_near = fit_params
        mises.SetParameters(amp_near, 0.0, kappa_near)
        mises.FixParameter(1, 0.0)

        dPhi.Fit(f'{self.name}_near_von_mises', 'BR')
    
        if ratio_plot:
            return rt.TRatioPlot(dPhi)
                
        else:
            return dPhi 
    
    
    def fit_away_von_mises(self, dPhi, fit_params, ratio_plot=False):
        '''
        Fit away-side with Von Mises. 
        '''
        dPhi.SetName(f'{self.name}_dPhi_fitted_away_mises')

        # 0: AS Amp/Yield, 1: AS mean, 2: AS kappa
        mises = rt.TF1(f"{self.name}_von_mises", "[0]/(2*TMath::Pi()*TMath::BesselI0([2]))*TMath::Exp([2]*TMath::Cos(x- 2*TMath::Pi() - [1]))",
                       rt.TMath.Pi()/2,
                       3*rt.TMath.Pi()/2)
        process_tf1(mises)
        mises.SetNpx(1000)

        amp_away, kappa_away = fit_params
        mises.SetParameters(amp_away, 0.0, kappa_away)
        mises.FixParameter(1, rt.TMath.Pi())

        dPhi.Fit(f'{self.name}_away_von_mises', 'B')
    
        if ratio_plot:
            return rt.TRatioPlot(dPhi)
                
        else:
            return dPhi 
        

    def fit_gen_gaussian(self, dPhi, fit_params, ratio_plot=False):
        '''
        Fit generalized Gaussian to delta phi. 
        '''
        
        dPhi.GetYaxis().SetRangeUser(0, 1.3 * dPhi.GetMaximum())
        #self._format_dist(dPhi)
        #process_hist(dPhi)
        dPhi.SetName(f'{self.name}_dPhi_fitted_gen')
        
        # 0: Mean NS, 1: alpha NS, 2: beta NS, 3: Amp GGaus NS, 4: Mean AS, 5: alpha AS, 6: beta AS, 7: Amp GGaus AS 
        genGauss = rt.TF1(f'{self.name}_fit', "[3]*([2]/(2*[1]*TMath::Gamma(1/[2])))*TMath::Exp(-TMath::Power(TMath::Abs(x-[0])/[1],[2])) + [7]*([6]/(2*[5]*TMath::Gamma(1/[6])))*TMath::Exp(-TMath::Power(TMath::Abs(x-[4])/[5],[6]))",
                          -rt.TMath.Pi()/2,
                          +3*rt.TMath.Pi()/2)
        process_tf1(genGauss)
        genGauss.SetNpx(1000)

        mean_near, alpha_near, beta_near, amp_near, mean_away, alpha_away, beta_away, amp_away = fit_params
        genGauss.SetParameters(mean_near, alpha_near, beta_near, amp_near, mean_away, alpha_away, beta_away, amp_away) 

        genGauss.FixParameter(0, 0.0)
        genGauss.FixParameter(4, rt.TMath.Pi())

        genGauss.SetParLimits(2, 1, 2)
        genGauss.SetParLimits(6, 1, 2)

        dPhi.Fit(f'{self.name}_fit', 'B')
        
        fit = dPhi.GetFunction(f'{self.name}_fit')
        params = np.array(list(zip([fit.GetParameter(i) for i in range(8)], [fit.GetParError(i) for i in range(8)])), dtype='d')  

        if ratio_plot:
            return rt.TRatioPlot(dPhi), params
                
        else:
            return dPhi, params
    
    
    @staticmethod
    def calc_gen_gaussian_widths(params):
        '''
        Calculate NS and AS standard deviations given fitted params from fit_gen_gaussian().
        '''
        alpha_near, alpha_near_uncert = params[1]
        beta_near, beta_near_uncert = params[2]

        alpha_away, alpha_away_uncert = params[5] 
        beta_away, beta_away_uncert = params[6]

        
        var_near = (alpha_near ** 2) * rt.TMath.Gamma(3 / beta_near) / rt.TMath.Gamma(1 / beta_near)
        var_away = (alpha_away ** 2) * rt.TMath.Gamma(3 / beta_away) / rt.TMath.Gamma(1 / beta_away)

        stdev_near = np.sqrt(var_near)
        stdev_away = np.sqrt(var_away)

        raise NotImplementedError('Implement uncertainties')
        return stdev_near, stdev_away


    def fit_near_gen_gaussian(self, dPhi, fit_params, ratio_plot=False):
        '''
        Fit generalized Gaussian to delta phi near-side. 
        '''
        
        dPhi.GetYaxis().SetRangeUser(0, 1.3 * dPhi.GetMaximum())
        #self._format_dist(dPhi)
        #process_hist(dPhi)
        dPhi.SetName(f'{self.name}_dPhi_near_fitted_gen')
        
        # 0: Mean NS, 1: alpha NS, 2: beta NS, 3: Amp GGaus NS
        genGauss = rt.TF1(f'{self.name}_fit', "[3]*([2]/(2*[1]*TMath::Gamma(1/[2])))*TMath::Exp(-TMath::Power(TMath::Abs(x-[0])/[1],[2]))",
                          -rt.TMath.Pi()/2,
                          rt.TMath.Pi()/2)
        process_tf1(genGauss)
        genGauss.SetNpx(1000)

        mean_near, alpha_near, beta_near, amp_near = fit_params
        genGauss.SetParameters(mean_near, alpha_near, beta_near, amp_near)
        genGauss.FixParameter(0, 0.0)
        genGauss.SetParLimits(2, 1, 2)

        dPhi.Fit(f'{self.name}_near_gen_gauss_fit', 'BR')

        if ratio_plot:
            return rt.TRatioPlot(dPhi)
                
        else:
            return dPhi
        
    
    def fit_away_gen_gaussian(self, dPhi, fit_params, ratio_plot=False):
        '''
        Fit generalized Gaussian to delta phi away-side. 
        '''
        
        dPhi.GetYaxis().SetRangeUser(0, 1.3 * dPhi.GetMaximum())
        #self._format_dist(dPhi)
        #process_hist(dPhi)
        dPhi.SetName(f'{self.name}_dPhi_away_fitted_gen')
        
        # 0: Mean AS, 1: alpha AS, 2: beta AS, 3: Amp GGaus AS 
        genGauss = rt.TF1(f'{self.name}_fit', "[3]*([2]/(2*[1]*TMath::Gamma(1/[2])))*TMath::Exp(-TMath::Power(TMath::Abs(x-[0])/[1],[2]))",
                          rt.TMath.Pi()/2,
                          3*rt.TMath.Pi()/2)
        process_tf1(genGauss)
        genGauss.SetNpx(1000)

        mean_away, alpha_away, beta_away, amp_away = fit_params
        genGauss.SetParameters(mean_away, alpha_away, beta_away, amp_away)
        genGauss.FixParameter(0, rt.TMath.Pi())
        genGauss.SetParLimits(2, 1, 2)

        dPhi.Fit(f'{self.name}_away_gen_gauss_fit', 'BR')

        if ratio_plot:
            return rt.TRatioPlot(dPhi)
                
        else:
            return dPhi
        
    def fit_near_side(self, dPhi, fit_params, ratio_plot=False):
        '''
        Fits the near-side peak of a delta phi distribution. The near side is defined as all points between a delta phi of -pi/2 and

        Args

        Returns

        '''

        dPhi.GetYaxis().SetRangeUser(0, 1.3 * dPhi.GetMaximum())
        #process_hist(dPhi)
        dPhi.SetName(f'{self.name}_dPhi_near_fitted')

        pi = rt.TMath.Pi()
        near_gaus = rt.TF1(f'{self.name}_near_fit', 'gaus', -pi/2, pi/2)
        process_tf1(near_gaus)
        near_gaus.SetNpx(1000)

        amp_near, stdev_near = fit_params
        near_gaus.SetParameters(amp_near, 0.0, stdev_near)
        near_gaus.FixParameter(1, 0.0)
        # near_gaus.SetParLimits(2, 0, 100)
        #near_gaus.FixParameter(2, 0.181)

        dPhi.Fit(f'{self.name}_near_fit', 'BR')

        
        fit = dPhi.GetFunction(f'{self.name}_near_fit')
        params = np.array(list(zip([fit.GetParameter(i) for i in range(3)], [fit.GetParError(i) for i in range(3)])), dtype='d')

        if ratio_plot:
            return rt.TRatioPlot(dPhi), params
                
        else:
            return dPhi, params
 

    def fit_away_side(self, dPhi, fit_params, ratio_plot=False):
        '''
        Args

        Returns
        '''
        dPhi.GetYaxis().SetRangeUser(0, 1.3 * dPhi.GetMaximum())
        #process_hist(dPhi)
        dPhi.SetName(f'{self.name}_dPhi_away_fitted')

        pi = rt.TMath.Pi()
        away_gaus = rt.TF1(f'{self.name}_away_fit', 'gaus', pi/2, 3*pi/2)
        process_tf1(away_gaus)
        away_gaus.SetNpx(1000)

        amp_away, stdev_away = fit_params
        away_gaus.SetParameters(amp_away, pi, stdev_away)
        away_gaus.FixParameter(4, pi)
        away_gaus.SetParLimits(5, 0, 100)

        dPhi.Fit(f'{self.name}_away_fit', 'BR')

        fit = dPhi.GetFunction(f'{self.name}_away_fit')
        params = np.array(list(zip([fit.GetParameter(i) for i in range(3)], [fit.GetParError(i) for i in range(3)])), dtype='d')
        #raise NotImplementedError
    
        if ratio_plot:
            return rt.TRatioPlot(dPhi), params
        else:
            return dPhi, params
        
        
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
        #self.fit_params = fit_params
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
            self.f = rt.TF1(f"{self.fit_name}", "[0]/(2*TMath::Pi()*TMath::BesselI0([2]))*TMath::Exp([2]*TMath::Cos(x- 2*TMath::Pi() - [1])) + [3]/(2*TMath::Pi()*TMath::BesselI0([5]))*TMath::Exp([5]*TMath::Cos(x- 2*TMath::Pi()-[4]))",
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