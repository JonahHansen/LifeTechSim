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


# =============================================================================
# SYSTEM
# =============================================================================

class System():
    
    def __init__(self,
                 Nuniverse,
                 Rp, # Rearth
                 Porb, # d
                 Mp, # Mearth
                 ep,
                 ip, # rad
                 Omegap, # rad
                 omegap, # rad
                 thetap, # rad
                 Abond,
                 AgeomVIS,
                 AgeomMIR,
                 z,
                 ap, # au
                 rp, # au
                 AngSep, # arcsec
                 maxAngSep, # arcsec
                 Fp, # Searth
                 fp,
                 Tp, # K
                 Nstar,
                 Rs, # Rsun
                 Ms, # Msun
                 Ts, # K
                 Ds, # pc
                 Stype,
                 RA, # deg
                 Dec): # deg
        """
        Parameters
        ----------
        Nuniverse: list
            Number of the universe to which the planet belongs to.
        Rp: list
            Planet radius (Rearth).
        Porb: list
            Planet orbital period (d).
        Mp: list
            Planet mass (Mearth).
        ep: list
            Planet eccentricity.
        ip: list
            Planet inclination (rad).
        Omegap: list
            Planet longitude of ascending node (rad).
        omegap: list
            Planet argument of periapsis (rad).
        thetap: list
            Planet true anomaly (rad).
        Abond: list
            Planet Bond albedo.
        AgeomVIS: list
            Planet geometric albedo in the visible.
        AgeomMIR: list
            Planet geometric albedo in the mid-infrared.
        z: list
            Exozodiacal dust level.
        ap: list
            Planet semi-major axis (au).
        rp: list
            Planet physical separation (au).
        AngSep: list
            Planet projected angular separation (arcsec).
        maxAngSep: list
            Max planet projected angular separation (arcsec).
        Fp: list
            Planet incident host star flux (Searth).
        fp: list
            Planet Lambertian reflectance.
        Tp: list
            Planet equilibrium temperature (K).
        Nstar: list
            Number of the star.
        Rs: list
            Host star radius (Rsun).
        Ms: list
            Host star mass (Msun).
        Ts: list
            Host star effective temperature (K).
        Ds: list
            Host star distance (pc).
        Stype: list
            Host star spectral type.
        RA: list
            Host star right ascension (deg).
        Dec: list
            Host star declination (deg).
        """
        
        self.Nuniverse = np.array(Nuniverse)
        self.Rp = np.array(Rp) # Rearth
        self.Porb = np.array(Porb) # d
        self.Mp = np.array(Mp) # Mearth
        self.ep = np.array(ep)
        self.ip = np.array(ip) # rad
        self.Omegap = np.array(Omegap) # rad
        self.omegap = np.array(omegap) # rad
        self.thetap = np.array(thetap) # rad
        self.Abond = np.array(Abond)
        self.AgeomVIS = np.array(AgeomVIS)
        self.AgeomMIR = np.array(AgeomMIR)
        self.z = np.array(z)
        self.ap = np.array(ap) # au
        self.rp = np.array(rp) # au
        self.AngSep = np.array(AngSep) # arcsec
        self.maxAngSep = np.array(maxAngSep) # arcsec
        self.Fp = np.array(Fp) # Searth
        self.fp = np.array(fp)
        self.Tp = np.array(Tp) # K
        self.Nstar = np.array(Nstar)
        self.Rs = np.array(Rs) # Rsun
        self.Ms = np.array(Ms) # Msun
        self.Ts = np.array(Ts) # K
        self.Ds = np.array(Ds) # pc
        self.Stype = np.array(Stype)
        self.RA = np.array(RA) # deg
        self.Dec = np.array(Dec) # deg
        
        pass
