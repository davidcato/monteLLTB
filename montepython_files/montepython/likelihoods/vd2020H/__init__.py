import os 
import numpy as np
from classy import Class
from matplotlib import *
import scipy.integrate as integrate
import scipy.constants as conts
import warnings
import montepython.io_mp as io_mp
import numexpr as ne

import scipy.linalg as la
from scipy.interpolate import interp1d
from sklearn.neighbors import KernelDensity
from scipy.optimize import minimize, fmin

from montepython.likelihood_class import Likelihood
from montepython.LLTB_functions import metric_functions, derivate_functions

class vd2020H(Likelihood):

    def __init__(self, path, data, command_line): # initialization routine

        Likelihood.__init__(self,path,data,command_line)

        ################# Importing SNe data
        dataSNE =  os.path.join(self.data_directory,self.lc_sn)
        covSNE =  os.path.join(self.data_directory,self.cov_sn)

        self.light_curve_params = self.read_light_curve_parameters_sn(dataSNE)
        self.C00  = self.read_matrix_sn(covSNE)

        dataSNE_low =  os.path.join(self.data_directory,self.lc_sn_low)
        covSNE_low =  os.path.join(self.data_directory,self.cov_sn_low)

        self.light_curve_params_low = self.read_light_curve_parameters_sn_low(dataSNE_low)
        self.C00_low  = self.read_matrix_sn(covSNE_low)

    def loglkl(self, cosmo, data, LLTBin):

        cc = conts.physical_constants['speed of light in vacuum'][0]/1000.

        ###################################
        # Cosmological and nuisance params
        ###################################
        MB = data.mcmc_parameters['MB']['current']

        VH0 = 100 * data.mcmc_parameters['h']['current'] * data.mcmc_parameters['h']['scale']
        Vomega_b = data.mcmc_parameters['omega_b']['current'] * data.mcmc_parameters['omega_b']['scale']
        Vomega_cdm = data.mcmc_parameters['omega_cdm']['current'] * data.mcmc_parameters['omega_cdm']['scale']
        VOmegaDE  = cosmo.Omega_Lambda()
        delta0 = data.mcmc_parameters['delta0']['current'] * data.mcmc_parameters['delta0']['scale']
        zB = data.mcmc_parameters['zB']['current'] * data.mcmc_parameters['zB']['scale']

        try:
            VOmegaK = data.mcmc_parameters['Omega_k']['current'] * data.mcmc_parameters['Omega_k']['scale']
        except:
            VOmegaK = 0.0

        # adding locals H0 as derived parameters
        data.derived_lkl = {'H0m':0.0,'H0l1':0.0,'H0l2':0.0,'H0R':0.0,'T0_eff':0.0,'L':0.0,'rB':0.0,'RL':0.0,'delta_L':0.0,'z_L':0.0}

        ###################################
        # Init LLTB and getting LLTB distances
        ###################################
        ini_distance = LLTBin
        LLTB_distance = ini_distance[1]
        nonsense_cosmo = ini_distance[0]
        params_eff = ini_distance[-1]
        LLTB_dic = {}

        ###################################
        # Computing likelihoods
        ###################################
        # skip likelihoods computation if nonsense_cosmo != 0
        # see likelihood_class.py
        if nonsense_cosmo == 0:

            ##################################
            # vd2020 distances and parameters
            ###################################
            # Dictionary with parameters used to compute
            # LLTB metric and derivate functions
            LLTB_dic['zB'] = zB
            LLTB_dic['VH0'] = VH0
            LLTB_dic['VOmegaDE'] = VOmegaDE
            LLTB_dic['VOmegaK'] = VOmegaK
            LLTB_dic['MB'] = MB

            mfun = metric_functions(LLTB_distance)
            dfun = derivate_functions(LLTB_distance,LLTB_dic)

            LoverB = params_eff['LoverB'][0]
            rBMpc = params_eff['rBMpc'][0]
            T0_eff = params_eff['Tcmb_eff'][0]

            ###################################
            # Local Hubble
            ###################################
            # Importing data
            redshifts_low = self.light_curve_params_low.zcmb
            dmb_low = self.light_curve_params_low.dmb
            C00_low = self.C00_low
            cov_low = ne.evaluate("C00_low")
            cov_low += np.diag(dmb_low**2)
            incov_low = la.inv(cov_low)

            # H0_m
            zsnlow = np.array(redshifts_low)
            zinthist = zsnlow[:, np.newaxis]
            zlen = len(zsnlow)
            zmin = np.min(zsnlow)
            zmax = np.max(zsnlow)
            zinter = np.linspace(zmin,zmax,5*zlen)

            kde = KernelDensity(kernel='epanechnikov', bandwidth=0.01).fit(zinthist)
            fhistz = kde.score_samples(zinter[:, np.newaxis])

            wz_sne_nnor = interp1d(zinter,np.exp(fhistz),fill_value='extrapolate')
            nnor = integrate.simps(np.exp(fhistz),zinter)

            wz_sne =interp1d(zinter,np.exp(fhistz)/nnor,fill_value='extrapolate')

            # H0mean --- eq 6 from WMC dLrel(self,z,rz,rb,dAz,OmDE,OmK,a0r,H0r)
            zinte = np.arange(zmin,zmax,0.0005)
            dLrelinte = dfun.dLrel(zinte)
            dLrelinte_weight = dLrelinte*wz_sne(zinte)
            H0m = cc/integrate.simps(dLrelinte_weight,zinte)

            # Mock LLTB supernovae
            mutab = dfun.muz_LLTB(zsnlow)

            # H0_l1
            # q0(r) renormalized by a constant hh
            def chi2_mu_l1(hh): 
                diff = dfun.mu_cosmo_LLTB(zsnlow,hh) - mutab
                chi2dL = np.dot(np.dot(diff,incov_low),diff)
                return chi2dL

            mymin_l1 = minimize(chi2_mu_l1, x0=[68.0])
            H0l1 = mymin_l1.x[0]

            # H0_l2
            # q0 and H0 free
            def chi2_mu_l2(vec): 
                hh,q0 = vec
                diff = dfun.mu_cosmo_FLRW(zsnlow,hh,q0) - mutab
                chi2dL = np.dot(np.dot(diff,incov_low),diff)
                return chi2dL

            mymin_l2 = minimize(chi2_mu_l2, x0=[68.0,-0.5], bounds=((10,120),(-4.,4.)))
            H0l2 = mymin_l2.x[0]

            # H0_R
            # H0 free and q0 = -0.55
            def chi2_mu_R(vec):
                hh = vec
                diff = dfun.mu_cosmo_FLRW(zsnlow,hh,-0.55) - mutab
                chi2dL = np.dot(np.dot(diff,incov_low),diff)
                return chi2dL

            mymin_R = minimize(chi2_mu_R, x0=[70.0])
            H0R = mymin_R.x[0]

            if zB == 0:
                zB = 1e-6

            rL = LoverB*mfun.rz_LLTB(zB)
            zL = mfun.zr_LLTB(rL)
            rLMpc = LoverB*rBMpc
            toMpc = rLMpc/rL
            RL = mfun.R0r_LLTB(rL)
            # RL is the physical distances
            deltaL = dfun.deltar_LLTB(RL)
            
            # Update derived parameter
            data.derived_lkl.update({'H0m':H0m,'H0l1':H0l1,'H0l2':H0l2,'H0R':H0R,'T0_eff':T0_eff,'L':LoverB,'rB':rBMpc,'RL':RL*toMpc,'delta_L':deltaL,'z_L':zL})
            
            # H0 likelihood using H0l1
            loglklH0 = -0.5 * (H0l1 - self.priorH0)**2 / (self.sigH0)**2

            ###################################
            # Likelihood SN
            ###################################
            # Importing data
            redshifts = self.light_curve_params.zcmb
            mb = self.light_curve_params.mb
            dmb = self.light_curve_params.dmb
            C00 = self.C00
            cov = ne.evaluate("C00")
            cov += np.diag(dmb**2)

            mb_LLTB = dfun.mB_LLTB(redshifts)
            difmB = mb - mb_LLTB 

            # Whiten the residuals, in two steps
            # 1) Compute the Cholesky decomposition of the covariance matrix, in
            # place. This is a time expensive (0.015 seconds) part
            cov = la.cholesky(cov, lower=True, overwrite_a=True)

            # 2) Solve the triangular system, also time expensive (0.02 seconds)
            difmB = la.solve_triangular(cov, difmB, lower=True, check_finite=False)

            # Finally, compute the chi2 as the sum of the squared residuals
            chi2sn = (difmB**2).sum()
            loglklsn = -0.5 * chi2sn

        else:
            # unphysical point
            data.derived_lkl.update({'H0m':-1,'H0l1':-1,'H0l2':-1,'H0R':-1,'T0_eff':-1,'L':-1.0,'rB':-1.0,'RL':-1.0,'delta_L':-1.0,'z_L':-1.0})
            loglklH0 = -1e12
            loglklsn = -1e12

        loglkl= loglklH0 + loglklsn # + another ones!

        return loglkl
