"""
# =============================================================================
# P-POP PHOTOMETRY
# A photometry tool for P-POP
# =============================================================================
"""


# =============================================================================
# IMPORTS
# =============================================================================

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simps


# =============================================================================
# BLACKBODY
# =============================================================================

class Photometry():
    
    def __init__(self):
        """
        """
        
        # Print.
        print('--> Initializing Blackbody')
        
        # Constants.
        self.h = 6.62607004e-34 # m^2*kg/s
        self.c = 299792458. # m/s
        self.kB = 1.38064852e-23 # m^2*kg/s^2/K
        self.Rsun = 695700000. # m
        self.pc = 3.0856776e16 # m
        
        pass
    
    def Compute(self,
                Filter,
                Sys,
                Unit,
                Mission):
        """
        Parameters
        ----------
        Filter: instance
            Instance of class Filter.
        Sys: instance
            Instance of class System.
        Unit: 'uJy', 'ph'
            Unit in which the photometry should be computed.
        Mission: 'MIR', 'VIS'
            Wavelength range in which the mission is operating.
        """
        
        if (Unit == 'uJy'):
            Flx = self.Flx_SI(Filter.Wavel, # m
                              Sys.Ts[0], # K
                              Sys.Rs[0], # Rsun
                              Sys.Ds[0]) # pc
            IntFlx = self.IntFlx_uJy(Flx, # W/m^3
                                     Filter.Wavel, # m
                                     Filter.Trans,
                                     Filter.Width, # m
                                     Filter.Mean) # m
        elif (Unit == 'ph'):
            Flx = self.Flx_ph(Filter.Wavel, # m
                              Sys.Ts[0], # K
                              Sys.Rs[0], # Rsun
                              Sys.Ds[0]) # pc
            IntFlx = self.IntFlx_ph(Flx, # ph/s/m^2/um
                                    Filter.Wavel, # m
                                    Filter.Trans,
                                    Filter.Width, # m
                                    Filter.Mean) # m
        
        return [IntFlx]*len(Sys.Nuniverse)
    
    def Flx_SI(self,
               Wavel, # m
               Ts, # K
               Rs, # Rsun
               Ds): # pc
        """
        Parameters
        ----------
        Wavel: array
            Wavelength (m) of filter nodes.
        Ts: float
            Host star effective temperature (K).
        Rs: float
            Host star radius (Rsun).
        Ds: float
            Host star distance (pc).
        
        Returns
        -------
        Flx: array
            Thermal blackbody flux (W/m^3).
        """
        
        Flx = 2.*np.pi*self.h*self.c**2/Wavel**5/(np.exp(self.h*self.c/(Wavel*self.kB*Ts))-1.)*((Rs*self.Rsun)/(Ds*self.pc))**2 # W/m^3
        
        return Flx
    
    def Flx_ph(self,
               Wavel, # m
               Ts, # K
               Rs, # Rsun
               Ds): # pc
        """
        Parameters
        ----------
        Wavel: array
            Wavelength (m) of filter nodes.
        Ts: float
            Host star effective temperature (K).
        Rs: float
            Host star radius (Rsun).
        Ds: float
            Host star distance (pc).
        
        Returns
        -------
        Flx: array
            Thermal blackbody flux (ph/s/m^3).
        """
        
        Flx = 2.*np.pi*self.c/Wavel**4/(np.exp(self.h*self.c/(Wavel*self.kB*Ts))-1.)*((Rs*self.Rsun)/(Ds*self.pc))**2 # ph/s/m^3
        
        return Flx
    
    def IntFlx_uJy(self,
                   Flx, # W/m^3
                   Wavel, # m
                   Trans,
                   Width, # m
                   Mean): # m
        """
        Parameters
        ----------
        Flx: array
            Thermal blackbody flux (W/m^3).
        Wavel: array
            Wavelength (m) of filter nodes.
        Trans: array
            Transmission of filter nodes.
        Width: float
            Width (m) of the filter.
        Mean: float
            Mean (m) of the filter.
        
        Returns
        -------
        Flx: array
            Integrated thermal blackbody flux (uJy).
        """
        
        IntFlx = 1e6*simps(Flx*Trans, Wavel)/Width*Mean**2/self.c*1e26 # uJy
        
        return IntFlx
    
    def IntFlx_ph(self,
                  Flx, # ph/s/m^3
                  Wavel, # m
                  Trans,
                  Width, # m
                  Mean): # m
        """
        Parameters
        ----------
        Flx: array
            Thermal blackbody flux (ph/s/m^3).
        Wavel: array
            Wavelength (m) of filter nodes.
        Trans: array
            Transmission of filter nodes.
        Width: float
            Width (m) of the filter.
        Mean: float
            Mean (m) of the filter.
        
        Returns
        -------
        Flx: array
            Integrated thermal blackbody flux (ph/s/m^2).
        """
        
        AbsTrans = 1.
        IntFlx = simps(Flx*Trans*AbsTrans, Wavel) # ph/s/m^2
        
        return IntFlx
    
    def SummaryPlots(self,
                     FigDir=None,
                     block=True):
        """
        Parameters
        ----------
        FigDir: str
            Directory to which summary plots are saved.
        block: bool
            If True, blocks plots when showing.
        """
        
        Wavel = np.logspace(-7, -5, 1000)
        Ts = 5800.
        Rs = 1.
        Ds = 149597870700./self.pc
        Flx_SI = self.Flx_SI(Wavel,
                             Ts,
                             Rs,
                             Ds)
        Flx_ph = self.Flx_ph(Wavel,
                             Ts,
                             Rs,
                             Ds)
        
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        plt.figure()
        ax0 = plt.gca()
        l0 = ax0.plot(Wavel*1e6, Flx_SI, color=colors[0], label='SI')
        ax0.set_xscale('log')
        ax0.set_yscale('log')
        ax0.grid(axis='y')
        ax0.set_xlabel('Wavelength [microns]')
        ax0.set_ylabel('Flux [W/m${}^3$]')
        ax1 = ax0.twinx()
        l1 = ax1.plot(Wavel*1e6, Flx_ph, color=colors[1], label='ph')
        ls = l0+l1
        la = [l.get_label() for l in ls]
        ax1.set_yscale('log')
        ax1.set_ylabel('Flux [ph/s/m${}^3$]', rotation=270, labelpad=20)
        ax1.legend(ls, la)
        plt.title('Sun @ 1 au')
        plt.tight_layout()
        if (FigDir is not None):
            plt.savefig(FigDir+'Blackbody.pdf')
        plt.show(block=block)
        plt.close()
        
        Wavel = np.linspace(1e-7, 1e-5, 1000000)
        Flx = self.Flx_SI(Wavel,
                          Ts,
                          Rs,
                          Ds)
        print('Solar irradiance @ 1 au = %.0f W/m^2' % simps(Flx, Wavel))
        
        pass
