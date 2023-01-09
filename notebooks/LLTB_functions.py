"""
.. module:: LLTB_functions
   :synopsis: Functions from LLTB model
.. moduleauthor:: David Camarena


Contains :class:`metric_functions` and :class:`derived_function`, with
basic functions, as well as more specific likelihood classes that may be reused
to implement new ones.

"""
import numpy as np
import scipy.constants as const
import scipy.interpolate as interpolate

class metric_functions(object):
    """
    Class with metric functions of LTB metric. 

    """

    def __init__(self, dataframe_LLTB):
        """
        Initialize the :class:`metric_functions`.

        Parameters
        ----------
        dataframe_LLTB : object
            `dataframe <pandas.core.frame.DataFrame>`
        """
        self.df = dataframe_LLTB

        try:
            self.df['z']
        except:
            raise ValueError('dataframe_LLTB seems to be empty.')

        self.rz_LLTB = interpolate.interp1d(self.df['z'],self.df['r'],fill_value='extrapolate')
        
        self.az_LLTB = interpolate.interp1d(self.df['z'],self.df['a'],fill_value='extrapolate')
        self.a0r_LLTB = interpolate.interp1d(self.df['r'],self.df['a0_r'],fill_value='extrapolate')
        
        self.dAz_LLTB = interpolate.interp1d(self.df['z'],self.df['d_A'],fill_value='extrapolate')
        
        self.Rz_LLTB = interpolate.interp1d(self.df['z'],self.df['R'],fill_value='extrapolate')
        self.Rztdrag_LLTB = interpolate.interp1d(self.df['z'],self.df['Rz_tdrag'],fill_value='extrapolate')
        self.HLz_LLTB = interpolate.interp1d(self.df['z'],self.df['HL'],fill_value='extrapolate')

        self.Rpztdrag_LLTB = interpolate.interp1d(self.df['z'],self.df['Rpz_tdrag'],fill_value='extrapolate') # In agreement with Rtdrag/r
        self.Rpdotz_LLTB = interpolate.interp1d(self.df['z'],self.df['Rpdot'],fill_value='extrapolate')

        self.H0r_LLTB = interpolate.interp1d(self.df['r'],self.df['H0_r'],fill_value='extrapolate')
        self.Hz_LLTB = interpolate.interp1d(self.df['z'],self.df['H'],fill_value='extrapolate')

    def Rpz_LLTB(self,z):
        rp = self.Rpdotz_LLTB(z)/self.HLz_LLTB(z)
        return rp

    def dLz_LLTB(self,z):
        dl = (1.0+z)**2*self.dAz_LLTB(z)
        return dl


class derivate_functions(metric_functions):
    """
    Class with derivate functions of LLTB model. 

    """

    def __init__(self, dataframe_LLTB, dict_LLTB_args,bao=False):
        """
        Initialize the :class:`metric_functions`.

        Parameters
        ----------
        dataframe_LLTB : object
            `dataframe <pandas.core.frame.DataFrame>`
        dict_LLTB_args : dictionary
            cosmological params needed to compute or evaluate derivate functions,
            zB, zdrag, rdrag, Omegas0 and MB are mandatory.
        """

        metric_functions.__init__(self, dataframe_LLTB)

        mandatory_params = ['zB','VH0','VOmegaDE','VOmegaK','MB'] #,'rs_eff','zdrag_eff'
        for manpa in mandatory_params:
            try:
                dict_LLTB_args[manpa]
            except:
                raise ValueError('%s is not included in dict_LLTB_args'%manpa)
        
        self.zboundary = dict_LLTB_args['zB']
        self.rboundary = self.rz_LLTB(self.zboundary)

        self.H = dict_LLTB_args['VH0']
        self.OmDE = dict_LLTB_args['VOmegaDE']
        self.OmK = dict_LLTB_args['VOmegaK']
        self.Mb = dict_LLTB_args['MB']

        if bao:
            self.rsd = dict_LLTB_args['rs_eff']
            self.zsd = dict_LLTB_args['zdrag_eff']
            self.lsd = self.rsd/(1.0+self.zsd )

    def muz_LLTB(self,z):
        mu = 5*np.log10(self.dLz_LLTB(z))+25
        return mu

    def mB_LLTB(self,z):
        mb=5*np.log10(self.dLz_LLTB(z)) +25 + self.Mb
        return mb

    def lT_drag(self,z):
        lT = self.Rz_LLTB(z)*self.lsd/self.Rztdrag_LLTB(z)
        return lT

    def lL_drag(self,z):
        lL = self.Rpz_LLTB(z)*self.lsd/self.Rpztdrag_LLTB(z)
        return lL

    def DeltaThetaz_LLTB(self,z):
        deltaThe = self.lT_drag(z)/self.dAz_LLTB(z)
        return deltaThe

    def DeltaZz_LLTB(self,z):
        deltazz = self.lL_drag(z)*(1.0+z)*self.HLz_LLTB(z)
        return deltazz

    def dVz_LLTB(self,z):
        cc = const.physical_constants['speed of light in vacuum'][0]/1000.
        invr = (self.DeltaZz_LLTB(z)*self.DeltaThetaz_LLTB(z)**2/(z*cc))**(1.0/3.0)
        dvz = self.rsd/invr
        return dvz

    def OmDEr_LLTB(self,r):
        omder = self.OmDE*self.H0r_LLTB(self.rboundary)**2/self.H0r_LLTB(r)**2
        return omder

    def OmMr_LLTB(self,r):
        ommr = (1.-self.OmDE-self.OmK)*self.a0r_LLTB(self.rboundary)**3/self.a0r_LLTB(r)**3*self.H0r_LLTB(self.rboundary)**2/self.H0r_LLTB(r)**2
        return ommr

    def OmKr_LLTB(self,r):
        omkr = 1.-self.OmDEr_LLTB(r)-self.OmMr_LLTB(r)
        return omkr

    def deltar_void(self,r):
        Omout = (1.-self.OmDE-self.OmK)
        dlt = -1 + self.OmMr_LLTB(r,rb,OmDE,OmK,a0r,H0r)*self.H0r_LLTB(r)**2/(Omout*H0r_LLTB(self.rboundary)**2)
        return dlt

    def q0r_LLTB(self,r):
        qzr = self.OmMr_LLTB(r)/2 - self.OmDEr_LLTB(r)
        return qzr

    def fcosmo_LLTB(self,z,hh):
        q0 = self.q0r_LLTB(self.rz_LLTB(z))*self.H0r_LLTB(self.rz_LLTB(z))**2/hh**2
        term0=1.0
        term1=(1.0-q0)/2.0
        fcos=term0*z+term1*z**2 
        return fcos

    def dL_cosmo_LLTB(self,z,hh):
        cc = const.physical_constants['speed of light in vacuum'][0]/1000.
        dl = cc*self.fcosmo_LLTB(z,hh)/hh
        return dl

    def mu_cosmo_LLTB(self,z,hh):
        mu = 5*np.log10(self.dL_cosmo_LLTB(z,hh))+25
        return mu

    def dLrel(self,z):
        drel = self.dLz_LLTB(z)/self.fcosmo_LLTB(z,self.H0r_LLTB(self.rz_LLTB(z)))
        return drel

    def fcosmo_FLRW(self,z,q0,j0=None,s0=None):
        term0=1.0
        term1=(1.0-q0)/2.0
        if j0 == None:
            term2 = 0.
        else:
            term2 = -(1.0-q0-3*q0**2+j0)/6.0
        if s0 == None:
            term3 = 0.
        else:
            term3 = (2-2*q0-15*q0**2+5*j0+10*q0*j0+s0)/24.0
        fcos=term0*z + term1*z**2 +term2*z**3 + term3*z**4
        return fcos

    def dL_cosmo_FLRW(self,z,hh,q0,j0=None,s0=None):
        cc = const.physical_constants['speed of light in vacuum'][0]/1000.
        dl = cc*self.fcosmo_FLRW(z,q0,j0,s0)/hh
        return dl

    def mu_cosmo_FLRW(self,z,hh,q0,j0=None,s0=None):
        mu = 5*np.log10(self.dL_cosmo_FLRW(z,hh,q0,j0,s0))+25
        return mu