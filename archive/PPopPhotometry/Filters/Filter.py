"""
# =============================================================================
# P-POP PHOTOMETRY
# A photometry tool for P-POP
# =============================================================================
"""


# =============================================================================
# IMPORTS
# =============================================================================

from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simps


# =============================================================================
# FILTER
# =============================================================================

class Filter():
    
    def __init__(self,
                 Name,
                 Wavel, # m
                 Trans,
                 Mean=None, # m
                 Width=None): # m
        """
        Returns
        -------
        Name: str
            Name of the filter.
        Wavel: array
            Wavelength (m) of filter nodes.
        Trans: array
            Transmission of filter nodes.
        Mean: float, None
            Mean (m) of the filter.
        Width: float, None
            Width (m) of the filter.
        """
        
        self.Name = Name
        self.Wavel = Wavel # m
        self.Trans = Trans
        self.Mean = Mean # m
        if (self.Mean is None):
            self.Mean = self.getMean() # m
        self.Width = Width # m
        if (self.Width is None):
            self.Width = self.getWidth() # m
        
        # Print.
        print('--> Initializing Filter '+self.Name)
        print('Mean wavelength = %.3f microns' % (self.Mean*1e6))
        print('Effective width = %.3f microns' % (self.Width*1e6))
        
        pass
    
    def getMean(self):
        """
        Returns
        -------
        Mean: float
            Mean (m) of the filter.
        """
        
        Mean = simps(self.Wavel*self.Trans, self.Wavel)/simps(self.Trans, self.Wavel) # m
        
        return Mean
    
    def getWidth(self):
        """
        Returns
        -------
        Width: float
            Width (m) of the filter.
        """
        
        Width = simps(self.Trans, self.Wavel)/np.max(self.Trans) # m
        
        return Width
    
    def SummaryPlot(self,
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
        
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        plt.figure()
        plt.plot(self.Wavel/1e-6, self.Trans, zorder=2)
        plt.axvline(self.Mean/1e-6, color=colors[1], label='Mean', zorder=1)
        Rect = Rectangle((self.Mean/1e-6-self.Width/2./1e-6, 0.), self.Width/1e-6, np.max(self.Trans), color=colors[2], alpha=1./3., label='Width', zorder=0)
        plt.gca().add_patch(Rect)
        plt.grid(axis='y')
        plt.xlabel('Wavelength [microns]')
        plt.ylabel('Transmission')
        plt.legend()
        plt.title(self.Name)
        plt.tight_layout()
        if (FigDir is not None):
            plt.savefig(FigDir+self.Name.replace('/', '_')+'.pdf')
        plt.show(block=block)
        plt.close()
        
        pass
