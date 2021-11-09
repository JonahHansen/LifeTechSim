"""
# =============================================================================
# P-POP PHOTOMETRY
# A photometry tool for P-POP
# =============================================================================
"""


# =============================================================================
# IMPORTS
# =============================================================================

import numpy as np
import sys

from . import System


# =============================================================================
# SYSTEMREADER
# =============================================================================

class SystemReader():

    def __init__(self,
                 PathPlanetTable):
        """
        Parameters
        ----------
        PathPlanetTable: str
            Path of the planet table to be read.
        """

        # Print.
        print('--> Initializing SystemReader')

        self.PathPlanetTable = PathPlanetTable

        pass

    def Open(self):
        """
        """

        # Open the planet table.
        Table = open(self.PathPlanetTable, 'r')
        self.Lines = Table.readlines()
        self.Nlines = len(self.Lines)
        Table.close()

        # The second line (i = 1) contains the column names of the new P-pop
        # while the first line (i = 0) contains the column names of the old
        # P-pop.
        tempLine = np.array(self.Lines[1].split('\t'))
        self.ColNuniverse = np.where(tempLine == 'Nuniverse')[0][0]
        self.ColRp = np.where(tempLine == 'Rp')[0][0]
        self.ColPorb = np.where(tempLine == 'Porb')[0][0]
        self.ColMp = np.where(tempLine == 'Mp')[0][0]
        self.Colep = np.where(tempLine == 'ep')[0][0]
        self.Colip = np.where(tempLine == 'ip')[0][0]
        self.ColOmegap = np.where(tempLine == 'Omegap')[0][0]
        self.Colomegap = np.where(tempLine == 'omegap')[0][0]
        self.Colthetap = np.where(tempLine == 'thetap')[0][0]
        self.ColAbond = np.where(tempLine == 'Abond')[0][0]
        self.ColAgeomVIS = np.where(tempLine == 'AgeomVIS')[0][0]
        self.ColAgeomMIR = np.where(tempLine == 'AgeomMIR')[0][0]
        self.Colz = np.where(tempLine == 'z')[0][0]
        self.Colap = np.where(tempLine == 'ap')[0][0]
        self.Colrp = np.where(tempLine == 'rp')[0][0]
        self.ColAngSep = np.where(tempLine == 'AngSep')[0][0]
        self.ColmaxAngSep = np.where(tempLine == 'maxAngSep')[0][0]
        self.ColFp = np.where(tempLine == 'Fp')[0][0]
        self.Colfp = np.where(tempLine == 'fp')[0][0]
        self.ColTp = np.where(tempLine == 'Tp')[0][0]
        self.ColNstar = np.where(tempLine == 'Nstar')[0][0]
        self.ColRs = np.where(tempLine == 'Rs')[0][0]
        self.ColMs = np.where(tempLine == 'Ms')[0][0]
        self.ColTs = np.where(tempLine == 'Ts')[0][0]
        self.ColDs = np.where(tempLine == 'Ds')[0][0]
        self.ColStype = np.where(tempLine == 'Stype')[0][0]
        self.ColRA = np.where(tempLine == 'RA')[0][0]
        self.ColDec = np.where(tempLine == 'Dec')[0][0]

        # Reset the line counter.
        self.Reset()

        pass

    def Reset(self):
        """
        """

        # The third line (i = 2) is the first line that contains planets.
        self.Counter = 2

        pass

    def Clear(self):
        """
        """

        self.Nuniverse = []
        self.Rp = [] # Rearth
        self.Porb = [] # d
        self.Mp = [] # Mearth
        self.ep = []
        self.ip = [] # rad
        self.Omegap = [] # rad
        self.omegap = [] # rad
        self.thetap = [] # rad
        self.Abond = []
        self.AgeomVIS = []
        self.AgeomMIR = []
        self.z = []
        self.ap = [] # au
        self.rp = [] # au
        self.AngSep = [] # arcsec
        self.maxAngSep = [] # arcsec
        self.Fp = [] # Searth
        self.fp = []
        self.Tp = [] # K
        self.Nstar = []
        self.Rs = [] # Rsun
        self.Ms = [] # Msun
        self.Ts = [] # K
        self.Ds = [] # pc
        self.Stype = []
        self.RA = [] # deg
        self.Dec = [] # deg

        pass

    def nextSystem(self):
        """
        Returns
        -------
        Sys: instance, None
            Instance of class System.
        """

        # Clear the system.
        self.Clear()

        # If there is no planet in the system yet or if the current planet
        # belongs to the same universe and the same star, add the current
        # planet to the system.
        if (self.Counter < self.Nlines):
            tempLine = self.Lines[self.Counter].split('\t')
        else:
            return None
        while (len(self.Nuniverse) == 0 or (self.Nuniverse[-1] == int(tempLine[self.ColNuniverse]) and self.Nstar[-1] == int(tempLine[self.ColNstar]))):
            self.Nuniverse += [int(tempLine[self.ColNuniverse])]
            self.Rp += [float(tempLine[self.ColRp])] # Rearth
            self.Porb += [float(tempLine[self.ColPorb])] # d
            self.Mp += [float(tempLine[self.ColMp])] # Mearth
            self.ep += [float(tempLine[self.Colep])]
            self.ip += [float(tempLine[self.Colip])] # rad
            self.Omegap += [float(tempLine[self.ColOmegap])] # rad
            self.omegap += [float(tempLine[self.Colomegap])] # rad
            self.thetap += [float(tempLine[self.Colthetap])] # rad
            self.Abond += [float(tempLine[self.ColAbond])]
            self.AgeomVIS += [float(tempLine[self.ColAgeomVIS])]
            self.AgeomMIR += [float(tempLine[self.ColAgeomMIR])]
            self.z += [float(tempLine[self.Colz])]
            self.ap += [float(tempLine[self.Colap])] # au
            self.rp += [float(tempLine[self.Colrp])] # au
            self.AngSep += [float(tempLine[self.ColAngSep])] # arcsec
            self.maxAngSep += [float(tempLine[self.ColmaxAngSep])] # arcsec
            self.Fp += [float(tempLine[self.ColFp])] # Searth
            self.fp += [float(tempLine[self.Colfp])]
            self.Tp += [float(tempLine[self.ColTp])] # K
            self.Nstar += [int(tempLine[self.ColNstar])]
            self.Rs += [float(tempLine[self.ColRs])] # Rsun
            self.Ms += [float(tempLine[self.ColMs])] # Msun
            self.Ts += [float(tempLine[self.ColTs])] # K
            self.Ds += [float(tempLine[self.ColDs])] # pc
            self.Stype += [str(tempLine[self.ColStype])]
            self.RA += [float(tempLine[self.ColRA])] # deg
            self.Dec += [float(tempLine[self.ColDec])] # deg
            self.Counter += 1
            if (self.Counter < self.Nlines):
                tempLine = self.Lines[self.Counter].split('\t')
            else:
                break

        # Create the system.
        Sys = System.System(self.Nuniverse,
                            self.Rp, # Rearth
                            self.Porb, # d
                            self.Mp, # Mearth
                            self.ep,
                            self.ip, # rad
                            self.Omegap, # rad
                            self.omegap, # rad
                            self.thetap, # rad
                            self.Abond,
                            self.AgeomVIS,
                            self.AgeomMIR,
                            self.z,
                            self.ap, # au
                            self.rp, # au
                            self.AngSep, # arcsec
                            self.maxAngSep, # arcsec
                            self.Fp, # Searth
                            self.fp,
                            self.Tp, # K
                            self.Nstar,
                            self.Rs, # Rsun
                            self.Ms, # Msun
                            self.Ts, # K
                            self.Ds, # pc
                            self.Stype,
                            self.RA, # deg
                            self.Dec) # deg

        sys.stdout.write('\r--> Planet %.0f of %.0f' % ((self.Counter-1), self.Nlines-2))
        sys.stdout.flush()

        return Sys
