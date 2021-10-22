import sys
sys.path.append("..")
import numpy as np
from PPop.ReadPlanetPopulation import PlanetPopulation
from engine.sim_functions import Planck_wrapper
from astropy import constants as const
from scipy.integrate import quad

sigma_sb = const.sigma_sb.value #Stefan-Boltzmann constant
h = const.h.value #Planck constant
c = const.c.value #Speed of light
k_B = const.k_B.value #Boltzmann constant
au = const.au.value #AU in m
pc = const.pc.value #parsec in m
R_sol = const.R_sun.value #Solar radius in m
R_earth = const.R_earth.value #Earth radii in m
L_sol = const.L_sun.value #Solar luminosity in W
rad2mas = np.degrees(1)*3600e3 #Number of milliarcsec in one radian

"""
Return luminosity in solar luminosity (from Stefan-Boltzmann law)
INPUTS:
    T = Temperature of star in K
    R = Radius of star in R_sol
OUTPUTS:
    Luminosity of star in L_sol
"""
def Luminosity(T,R):
    Lum = 4*np.pi*(R*R_sol)**2*sigma_sb*T**4
    return Lum/L_sol


"""
Classes for Stars and Planets, to change the "ReadPlanetPopulation" data format
into something more usable
"""
class Star():
    def __init__(self,
                 Name,
                 Number,
                 Dist, # pc
                 Stype,
                 Rad, # Rsun
                 Teff, # K
                 Mass, # Msun
                 RA, # deg
                 Dec, # deg
                 spectrograph, #Class holding spectral info
                 HzMin, #au
                 HzMax, #au
                 z #exozodis
                 ):

        self.Name = Name #Name of star
        self.SNumber = Number #Number(index) of star
        self.Stype = Stype #Stellar type

        self.Dist = Dist #Distance in pc
        self.SRad = Rad #Radius in Rsol
        self.STeff = Teff #Effective Temperature in K
        self.SMass = Mass #Mass in Msol

        self.RA = RA #Ra in deg
        self.Dec = Dec #Dec in deg

        #Function for the stellar flux as a function of wavelength
        flux_function = self.star_flux_func()

        #Integrate over channels
        flux_ls = []
        for wave in spectrograph.channel_borders:
            flux_ls.append(quad(flux_function,wave,wave+spectrograph.dlambda)[0])
        self.flux = np.array(flux_ls)

        self.HZMin = HzMin #Outer? limit on HZ
        self.HZMax = HzMax #Inner? limit on HZ
        self.HZAngle = (0.5*(HzMin + HzMax)*au/(Dist*pc))*rad2mas #Average projected angle of HZ

        self.Exzod = z #Exozodiacal light level

        self.SLum = Luminosity(Teff,Rad) #Luminosity in Lsol

        parallax_mas = 1000/Dist

        self.HZEst = np.sqrt(self.SLum)*parallax_mas #Estimate of HZ from Luminosity

        self.angRad = Rad/215*parallax_mas #angular radius of star in mas

        self.Planets = [] #List of planets associated with the star (multiple universes)

    def star_flux_func(self):
        #Function to calculate the flux of star at a given wavelength
        R = self.SRad*R_sol #From solar radii to m
        d = self.Dist*pc #from parsec to meters
        p = Planck_wrapper(self.STeff)
        const = (R/d)**2
        def func(lam): #lam in m
            B = const*p(lam) #photons/m^2/s/m
            return B
        return func


class Planet():
    def __init__(self,
                 star,
                 UNumber,
                 SNumber,
                 PNumber,
                 PRad, #REarth
                 PMass, #MEarth
                 Porb, #days
                 e,
                 i, #rad
                 Omega, #rad
                 omega, #rad
                 theta, #rad
                 Abond,
                 AgeomVIS,
                 AgeomMIR,
                 a, #au
                 rSep, #au
                 AngSep, #arcsec
                 maxAngSep, #arcsec
                 F, #SEarth
                 f,
                 T, #K
                 spectrograph
                 ):

        self.parent_star = star

        self.Name = "U%dS%dP%d"%(UNumber,SNumber,PNumber)
        self.UNumber = UNumber #Universe index
        self.PNumber = PNumber #Planet index (per star, universe)

        self.PRad = PRad #Planet radius in REarth
        self.PMass = PMass #Planet mass in MEarth
        self.POrb = Porb #Planet orbital period in days

        self.e = e #Eccentricity
        self.i = i #Inclination in deg
        self.Omega = Omega #Longitude of ascending node in deg
        self.omega = omega #Argument of periapsis in deg
        self.theta = theta #True anomaly in deg

        self.Abond = Abond #Bond albedo
        self.AgeomVIS = AgeomVIS #Geometric visible albedo
        self.AgeomMIR = AgeomMIR #Geometric MIR albedo

        self.a = a #Semi major axis in au
        self.PrSep = rSep #Physical separation of planet in au
        self.PAngSep = AngSep*1000 #Projected angular sep in milliarcsec
        self.PmaxAngSep = maxAngSep*1000 #Max projected angular sep in milliarcsec

        self.IncFlux = F #Incident stellar flux in SEarth
        self.LamRef = f #Lambertian reflectance
        self.PTemp = T #Planet equilibrium temperature in K

        #Functions to get total planet flux as a function of wavelength
        reflected = self.reflected_flux_func()
        thermal = self.thermal_flux_func()
        total_flux = lambda lam: reflected(lam) + thermal(lam)

        #Integrate over wavelength
        flux_ls = []
        for wave in spectrograph.channel_borders:
            flux_ls.append(quad(total_flux,wave,wave+spectrograph.dlambda)[0])
        self.flux = np.array(flux_ls)

    #SEE PPOPPHOTOMETRY FOR THESE TWO FUNCTIONS
    def reflected_flux_func(self):
        #Function to calculate the flux of star at a given wavelength
        Rs = self.parent_star.SRad*R_sol #From solar radii to m
        d = self.parent_star.Dist*pc #from parsec to meters
        Rp = self.PRad*R_earth #From solar radii to m
        a = self.a*au #From au to m
        p = Planck_wrapper(self.parent_star.STeff)
        const = self.AgeomMIR*self.LamRef*(Rs/d)**2*(Rp/a)**2
        def func(lam): #lam in m
            B = const*p(lam) #photons/m^2/s/m
            return B
        return func

    def thermal_flux_func(self):
        #Function to calculate the thermal flux of planet at a given wavelength
        R = self.PRad*R_earth #From solar radii to m
        d = self.parent_star.Dist*pc #from parsec to meters
        p = Planck_wrapper(self.PTemp)
        const = (R/d)**2
        def func(lam): #lam in m
            B = const*p(lam) #photons/m^2/s/m
            return B
        return func

"""
Take star and planet data from PPop and return a list of star/planet datastructures
to use in simulation
INPUTS:
    planet_path = path to planet data
    phot_path = path to photometry data
OUTPUT:
    List of star datastructures, that each contain a list of planets
    ([planets] = Star[i].Planets)

"""
def RetrievePlanetData(planet_path,spectrograph):

    #Get planet data from PPop
    planet_data = PlanetPopulation(planet_path)
    #Get HZ data
    planet_data.ComputeHZ(Model='MS')
    #Attach photometry
    #planet_data.appendPhotometry(phot_path,"tag")

    star_list = [] #List of stars
    planet_list = [] #List of planets (for each star)
    universe_no = planet_data.Nuniverse[0] #Universe index
    star_no = planet_data.Nstar[0] #Star index
    planet_no = 0 #Planet index

    #Initial star
    star = Star(planet_data.SName[0],
                planet_data.Nstar[0],
                planet_data.Ds[0],
                planet_data.Stype[0],
                planet_data.Rs[0],
                planet_data.Ts[0],
                planet_data.Ms[0],
                planet_data.RA[0],
                planet_data.Dec[0],
                spectrograph,
                planet_data.HZin[0],
                planet_data.HZout[0],
                planet_data.z[0]
                )

    #Add star to list
    star_list.append(star)

    #Go through data structure
    for i in range(len(planet_data.Rp)):

        if planet_data.Nstar[i] != star_no: #New star
            #Attach planets to previous star
            star_list[-1].Planets = planet_list
            planet_list = [] #Reset list
            star_no = planet_data.Nstar[i] #Next index

            #New star data
            star = Star(planet_data.SName[i],
                        planet_data.Nstar[i],
                        planet_data.Ds[i],
                        planet_data.Stype[i],
                        planet_data.Rs[i],
                        planet_data.Ts[i],
                        planet_data.Ms[i],
                        planet_data.RA[i],
                        planet_data.Dec[i],
                        spectrograph,
                        planet_data.HZin[i],
                        planet_data.HZout[i],
                        planet_data.z[i]
                        )
            star_list.append(star)

        #If we are in the next universe
        if planet_data.Nuniverse[i] != universe_no:
            universe_no = planet_data.Nuniverse[i] #Next index
            planet_no = 0 #Reset planet index

        #Planet data
        planet = Planet(star,
                        planet_data.Nuniverse[i],
                        planet_data.Nstar[i],
                        planet_no,
                        planet_data.Rp[i],
                        planet_data.Mp[i],
                        planet_data.Porb[i],
                        planet_data.ep[i],
                        planet_data.ip[i],
                        planet_data.Omegap[i],
                        planet_data.omegap[i],
                        planet_data.thetap[i],
                        planet_data.Abond[i],
                        planet_data.AgeomVIS[i],
                        planet_data.AgeomMIR[i],
                        planet_data.ap[i],
                        planet_data.rp[i],
                        planet_data.AngSep[i],
                        planet_data.maxAngSep[i],
                        planet_data.Fp[i],
                        planet_data.fp[i],
                        planet_data.Tp[i],
                        spectrograph)

        planet_list.append(planet)
        planet_no += 1

    #Attach planets to final star
    star_list[-1].Planets = planet_list

    return star_list
