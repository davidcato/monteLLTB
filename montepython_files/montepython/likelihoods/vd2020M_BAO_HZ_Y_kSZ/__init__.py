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

class vd2020M_BAO_HZ_Y_kSZ(Likelihood):

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

        ################ Importing BAO data (BAO-lowz and BOSS-DR12 consensus)
        # BOSS-DR12
        # define arrays for values of z and data points
        self.z_DR12 = np.array([], 'float64')
        self.DM_rdfid_by_rd_in_Mpc = np.array([], 'float64')
        self.H_rd_by_rdfid_in_km_per_s_per_Mpc = np.array([], 'float64')

        # read redshifts and data points
        with open(os.path.join(self.data_directory, self.data_file_DR12), 'r') as filein:
            for i, line in enumerate(filein):
                if line.strip() and line.find('#') == -1:
                    this_line = line.split()
                    # load redshifts and D_M * (r_s / r_s_fid)^-1 in Mpc
                    if this_line[1] == 'dM(rsfid/rs)':
                        self.z_DR12 = np.append(self.z_DR12, float(this_line[0]))
                        self.DM_rdfid_by_rd_in_Mpc = np.append(
                            self.DM_rdfid_by_rd_in_Mpc, float(this_line[2]))
                    # load H(z) * (r_s / r_s_fid) in km s^-1 Mpc^-1
                    elif this_line[1] == 'Hz(rs/rsfid)':
                        self.H_rd_by_rdfid_in_km_per_s_per_Mpc = np.append(
                            self.H_rd_by_rdfid_in_km_per_s_per_Mpc, float(this_line[2]))

        # read covariance matrix
        self.cov_data_DR12 = np.loadtxt(os.path.join(self.data_directory, self.cov_file_DR12))
        # number of bins
        self.num_bins = np.shape(self.z_DR12)[0]

        # BAO-lowz        
        # define array for values of z and data points
        self.z_lowz = np.array([], 'float64')
        self.data_lowz = np.array([], 'float64')
        self.error_lowz = np.array([], 'float64')
        self.type_lowz = np.array([], 'int')

        exclude_lowz = []

        # read redshifts and data points
        with open(os.path.join(self.data_directory, self.data_file_lowz), 'r') as filein:
            for line in filein:
                if line.strip() and line.find('#') == -1:
                    # the first entry of the line is the identifier
                    this_line = line.split()
                    # insert into array if this id is not manually excluded
                    if not this_line[0] in exclude_lowz:
                        self.z_lowz = np.append(self.z_lowz, float(this_line[1]))
                        self.data_lowz = np.append(self.data_lowz, float(this_line[2]))
                        self.error_lowz = np.append(self.error_lowz, float(this_line[3]))
                        self.type_lowz = np.append(self.type_lowz, int(this_line[4]))

        # number of data points
        self.num_points_lowz = np.shape(self.z_lowz)[0]

        ################ Importing Hz data
        # Read the content of the data file, containing z, Hz and error
        total = np.loadtxt(os.path.join(self.data_directory_hz, self.data_file_hz))

        # Store the columns separately
        self.z_hz = total[:, 0]
        self.Hz = total[:, 1]
        self.err_hz = total[:, 2]

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

        # needed for y-Compton and kSZ effect
        zreio = cosmo.get_current_derived_parameters(['z_reio'])['z_reio']
        YHe = cosmo.get_current_derived_parameters(['YHe'])['YHe']

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

        # effective FLRW to compute rs_drag
        if nonsense_cosmo == 0:

            h_eff = params_eff['H0_eff'][0]/100.0
            omega_b_eff = params_eff['omegah2_b_eff'][0]
            omega_dm_eff = params_eff['omegah2_dm_eff'][0]
            OmegaK_eff = params_eff['OmegaK_eff'][0]
            T0_eff = params_eff['Tcmb_eff'][0]

            n_s_in = data.mcmc_parameters['n_s']['current'] * data.mcmc_parameters['n_s']['scale']
            logA_s_in = data.mcmc_parameters['ln10^{10}A_s']['current'] * data.mcmc_parameters['ln10^{10}A_s']['scale']
            tau_reio_in = data.mcmc_parameters['tau_reio']['current'] * data.mcmc_parameters['tau_reio']['scale']
            N_ur=cosmo.Neff()

            newcosmo_arguments = {'h': h_eff, 'Omega_k': OmegaK_eff, 'tau_reio': tau_reio_in, 'omega_cdm': omega_dm_eff, 'omega_b': omega_b_eff, 'T_cmb': T0_eff, 
                                  'N_ur': N_ur, 'ln10^{10}A_s': logA_s_in, 'n_s': n_s_in}
    	
    	    from classy import Class
            newcosmo = Class()
            newcosmo.set(newcosmo_arguments)

            newcosmo.compute()
            rs_eff = newcosmo.rs_drag()
            zdrag_eff = newcosmo.get_current_derived_parameters(['z_d'])['z_d']
            LLTB_dic['rs_eff'] = rs_eff
            LLTB_dic['zdrag_eff'] = zdrag_eff


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
            dfun = derivate_functions(LLTB_distance,LLTB_dic,bao=True)

            LoverB = params_eff['LoverB'][0]
            rBMpc = params_eff['rBMpc'][0]

            # Concerning betaz and its derivates we keep the minimalist approach. Instead of interpolated functions
            # we simply use arrays. Remember that the goal is an integral.
            LLTB_distance['ycompton_inte'] = LLTB_distance['betaz']**2*LLTB_distance['dtaudrz']
            # To ensure that we are integrating from zero to zreio
            LLTB_distance['ycompton_inte'][LLTB_distance['z'] > zreio]=0

            ###################################
            # MB prior
            ###################################
            loglklMB = -0.5 * (MB - self.priorMB)**2 / (self.sigMB)**2

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

            ###################################
            # Likelihood BAO
            ###################################
            # BOSS-DR12
            # define array for  values of D_M_diff = D_M^th - D_M^obs and H_diff = H^th - H^obs,
            # ordered by redshift bin (z=[0.38, 0.51, 0.61]) as following:
            # data_array = [DM_diff(z=0.38), H_diff(z=0.38), DM_diff(z=0.51), .., .., ..]
            # zdrag = cosmo.get_current_derived_parameters(['z_d'])['z_d']
            # rd = cosmo.rs_drag() * self.rs_rescale_DR12

            data_array = np.array([], 'float64')
            # for each point, compute comoving angular diameter distance D_M = (1 + z) * D_A,
            # sound horizon at baryon drag rs_d, theoretical prediction
            for i in range(self.num_bins):
                # replace with LLTB 
                # DM_at_z = cosmo.angular_distance(self.z_DR12[i]) * (1. + self.z_DR12[i])
                # H_at_z = cosmo.Hubble(self.z_DR12[i]) * conts.c / 1000.0
                del_theta_LLTB_at_z = dfun.DeltaThetaz_LLTB(self.z_DR12[i])
                del_z_LLTB_at_z = dfun.DeltaZz_LLTB(self.z_DR12[i])

                # theo_DM_rdfid_by_rd_in_Mpc = DM_at_z / rd * rd_fid_in_Mpc_DR12 
                # theo_H_rd_by_rdfid = H_at_z * rd / rd_fid_in_Mpc_DR12 
                theo_DM_rdfid_by_rd_in_Mpc =  self.rd_fid_in_Mpc_DR12 / del_theta_LLTB_at_z
                theo_H_rd_by_rdfid = del_z_LLTB_at_z / self.rd_fid_in_Mpc_DR12 

                # calculate difference between the sampled point and observations
                DM_diff = theo_DM_rdfid_by_rd_in_Mpc - self.DM_rdfid_by_rd_in_Mpc[i]
                H_diff = theo_H_rd_by_rdfid - self.H_rd_by_rdfid_in_km_per_s_per_Mpc[i]

                # save to data array
                data_array = np.append(data_array, DM_diff)
                data_array = np.append(data_array, H_diff)

            # compute chi squared
            inv_cov_data_DR12 = np.linalg.inv(self.cov_data_DR12)
            chi2DR12 = np.dot(np.dot(data_array,inv_cov_data_DR12),data_array)

            # return ln(L)
            loglklDR12 = - 0.5 * chi2DR12

            # BAO-lowz
            chi2lowz = 0.

            # for each point, compute angular distance da, radial distance dr,
            # volume distance dv, sound horizon at baryon drag rs_d,
            # theoretical prediction and chi2 contribution
            for i in range(self.num_points_lowz):

                dv = dfun.dVz_LLTB(self.z_lowz[i])
                # rs should be computed using effective params
                # rs = cosmo.rs_drag()
                rs = rs_eff

                if self.type_lowz[i] == 3: #Use this
                    theo = dv / rs

                elif self.type_lowz[i] == 4:
                    theo = dv

                elif self.type_lowz[i] == 5:
                    theo = da / rs

                elif self.type_lowz[i] == 6: #Here is not possible
                    theo = 1. / cosmo.Hubble(self.z_lowz[i]) / rs

                elif self.type_lowz[i] == 7: #Use this too
                    theo = rs / dv
                else:
                    raise io_mp.LikelihoodError(
                        "In likelihood %s. " % self.name +
                        "BAO data type %s " % self.type_lowz[i] +
                        "in %d-th line not understood" % i)

                chi2lowz += ((theo - self.data_lowz[i]) / self.error_lowz[i]) ** 2

            # return ln(L)
            lkllowz = - 0.5 * chi2lowz

            loglklBAO = loglklDR12 + lkllowz

            # cleaning up effective FLRW
            newcosmo.struct_cleanup()
            newcosmo.empty()

            ###################################
            # Cosmic clocks Hz
            ###################################
            # Loop over the redshifts
            chi2hz = 0.
            for index, z in enumerate(self.z_hz):
                # Query the cosmo module for the Hubble rate (in 1/Mpc), and
                # convert it to km/s/Mpc
                H_cosmo = mfun.HLz_LLTB(z)
                # Add to the tota chi2
                chi2hz += (self.Hz[index]-H_cosmo)**2/self.err_hz[index]**2

            loglklhz = -0.5 * chi2hz

            ###################################
            # Likelihood Y Compton
            ###################################
            # Trapezoidal and Simpson integral.
            ytrap = 0.7*integrate.trapz(LLTB_distance['ycompton_inte'],LLTB_distance['z'])
            ysimp = 0.7*integrate.simps(LLTB_distance['ycompton_inte'],LLTB_distance['z'])
            dely = (1-(ysimp+1e-8)/(ytrap+1e-81)) # 1e-8 offset avoids divergences

            # Sanity check. If the diferences between both integral is greater than 10% computation
            # can be consider inconsistent
            if dely > 0.1:
                warnings.warn('Something wrong with Y Compton integral: '+r'$\delta_0 = $'+'%s'%delta0+' and '+r'$z_B = $'+'%s'%zB)

            if ysimp < self.upY:
                loglklY = 0.0
            else:
                loglklY = -1e12

            ###################################
            # Likelihood kSZ
            ###################################
            # we integrate from r_min to r_reio, where r_min correspond to kmax scale. 
            # We impose kmax <= 200
            kmax_try = np.max(LLTB_distance['k_back'])
            kmax_cut = 200.

            if kmax_cut < kmax_try:
                # cutting scales bigger than kmax_cut
                k_pass = np.array(LLTB_distance['k_back'][LLTB_distance['k_back'] <= kmax_cut])
                z_pass = np.array(LLTB_distance['z'][LLTB_distance['k_back'] <= kmax_cut])
                r_pass = np.array(LLTB_distance['r'][LLTB_distance['k_back'] <= kmax_cut])
                LLTB_distance['prefac_inte'] = LLTB_distance['r']*(LLTB_distance['betaz']*LLTB_distance['dtaudrz'])**2
                prefac_inte = np.array(LLTB_distance['prefac_inte'][LLTB_distance['k_back'] <= kmax_cut])
                # integral goes up to zreio
                k_pass = k_pass[z_pass<=zreio]
                r_pass = r_pass[z_pass<=zreio]
                prefac_inte = prefac_inte[z_pass<=zreio]
                z_pass = z_pass[z_pass<=zreio]
            else:
                # kmax_try scale is smaller than kmax_cut
                # integral goes up to zreio
                k_pass = np.array(LLTB_distance['k_back'][LLTB_distance['z']<=zreio][1:])
                z_pass = np.array(LLTB_distance['z'][LLTB_distance['z']<=zreio][1:])
                r_pass = np.array(LLTB_distance['r'][LLTB_distance['z']<=zreio][1:])
                LLTB_distance['prefac_inte'] = LLTB_distance['r']*(LLTB_distance['betaz']*LLTB_distance['dtaudrz'])**2
                prefac_inte = np.array(LLTB_distance['prefac_inte'][LLTB_distance['z']<=zreio][1:])

            if rBMpc != 0: # is not LCDM
                # avoid error on Class when there exists blueshift (for instance  delta0 = 0.9718786882631633 & zB = 0.08849359160146264)
                rB = mfun.rz_LLTB(zB)
                factorMpc = rBMpc/rB
                # use z_LCDM instead of z_LLTB
                zr_LCDM = interp1d(cosmo.get_background()[u'comov. dist.']/factorMpc,cosmo.get_background()[u'z'],kind='cubic')
                new_z_pass = zr_LCDM(r_pass)
                z_pass = new_z_pass
                              
            # getting kmax
            kmax = np.max(k_pass)

            if np.min(z_pass[0]) < 0:
                warnings.warn('Something wrong with kSZ, negative z : '+r'$\delta_0 = $'+'%s'%delta0+' and '+r'$z_B = $'+'%s'%zB)
                z_pass = np.abs(z_pass)

            # Same background cosmology with Power spectrum. Note that kmax is defined after call vd2020
            newcosmo_arguments = {'h': VH0/100., 'Omega_k': 0., 'tau_reio': tau_reio_in, 'omega_cdm': Vomega_cdm, 'omega_b': Vomega_b,
                                  'N_ur': N_ur, 'ln10^{10}A_s': logA_s_in, 'n_s': n_s_in, 'output': 'mPk','non linear':'Halofit',
                                  'z_max_pk': 15.0, 'P_k_max_h/Mpc': 1.5*kmax}

            newcosmo_arguments_high = {'h': VH0/100., 'Omega_k': 0., 'tau_reio': tau_reio_in, 'omega_cdm': Vomega_cdm, 'omega_b': Vomega_b,
                                  'N_ur': N_ur, 'ln10^{10}A_s': logA_s_in, 'n_s': n_s_in, 'output': 'mPk','non linear':'Halofit',
                                  'z_max_pk': 15.0, 'P_k_max_h/Mpc': 3.75*kmax}

            newcosmo = Class()
            newcosmo.set(newcosmo_arguments)

            newcosmo.compute()
            pk_prove = np.vectorize(newcosmo.pk)

            try:
                pk_prove(k_pass,z_pass)
            except:
                warnings.warn('Increasing k_max')
                newcosmo.struct_cleanup()
                newcosmo.empty()
                newcosmo.set(newcosmo_arguments_high)
                newcosmo.compute()

            sig8 = newcosmo.get_current_derived_parameters(['sigma8'])['sigma8']

            def pk_dimless(k,z):
                pk = k**3*newcosmo.pk(k,z)/(2*np.pi)
                return pk

            pk_dimless_back = np.vectorize(pk_dimless)
            delPk_sq_inte = pk_dimless_back(k_pass,z_pass)
            inte_tot = prefac_inte*delPk_sq_inte

            # eq. from VMC
            ClkSZ = 16.* integrate.trapz(inte_tot,r_pass) * np.pi**2 * T0_eff**2 * 1e12/(2*3000+1)**3
            DkSZ = ClkSZ * 3000 * 3001/(2*np.pi)

            # reio parametrizate as footnote 14 of 1807.06209 or eq. B3 of 0804.3865
            # so Dzreio = z25p - z75p can be computed once we know zreio
            # {z25p, z75p} = (1+zreio)[3/2*wireio/(1+zreio)*arctanh({-1,1}/2)-1]^(1/expreio) - 1
            # where wireio = reionization_width = 0.5 and expreio = reionization_exponent = 1.5 
            # expreio = 1.5
            wireio = 0.5
            z25p = (1.+zreio)*((1.5*wireio*np.arctanh(-0.5)/(1.+zreio) - 1)**2.)**(1./3.) - 1
            z75p = (1.+zreio)*((1.5*wireio*np.arctanh(0.5)/(1.+zreio) - 1)**2.)**(1./3.) - 1
            Dzreio = z25p - z75p
            # eqs 5 and 6 from 1406.4794
            hDkSZ = 1.65*(sig8/0.8)**4.46
            pDkSZ = 2.03*((1.+zreio)/11.-0.12)*(Dzreio/1.05)**0.51
            totDkSZ = DkSZ + hDkSZ + pDkSZ

            newcosmo.struct_cleanup()
            newcosmo.empty()

            loglklkSZ = -0.5 * (totDkSZ - self.priorkSZ)**2 / (self.sigkSZ)**2

        else:
            # unphysical point
            data.derived_lkl.update({'H0m':-1,'H0l1':-1,'H0l2':-1,'H0R':-1,'T0_eff':-1,'L':-1.0,'rB':-1.0,'RL':-1.0,'delta_L':-1.0,'z_L':-1.0})
            loglklMB = -1e12
            loglklsn = -1e12 
            loglklBAO = -1e12
            loglklhz = -1e12
            loglklY = -1e12
            loglklkSZ = -1e12

        loglkl= loglklMB + loglklsn + loglklBAO + loglklhz + loglklY + loglklkSZ # + another ones!

        return loglkl
