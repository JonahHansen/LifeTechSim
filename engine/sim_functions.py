import sys
sys.path.append("..")
import numpy as np
from jwst_backgrounds import jbt
from scipy.interpolate import interp1d
from astropy import constants as const

h = const.h.value #Planck constant
c = const.c.value #Speed of light
k_B = const.k_B.value #Boltzmann constant
au = const.au.value #AU in m
pc = const.pc.value #parsec in m
rad2mas = np.degrees(1)*3600e3 #Number of milliarcsec in one radian

#Initiate units in um. Attributes in m
class Spectrograph():
    def __init__(self, wave_min, wave_max, baseline_wave, n_channels):

        #Bandwidth of spectrograph
        self.bandwidth = (wave_max-wave_min)*1e-6
        #Mean wavelength of spectrograph
        self.mean = (wave_max+wave_min)/2*1e-6
        #Min and max wavelengths
        self.wave_min = wave_min*1e-6
        self.wave_max = wave_max*1e-6
        #Number of channels
        self.n_channels = n_channels

        #Wavelength to set the baseline (optimisation baseline)
        self.baseline_wave = baseline_wave*1e-6

        #Channel borders (start of each channel)
        self.channel_borders = np.linspace(self.wave_min,self.wave_max,n_channels+1)[:-1]
        #Size of channel
        self.dlambda = self.channel_borders[1]-self.channel_borders[0]
        #Channel centres (middle of each channel)
        self.channel_centres = (self.channel_borders + self.dlambda/2)
        #Effective resolution of the spectrograph
        self.eff_resolution = self.mean/self.dlambda




"""
Planck function
Return a function for the spectral flux density as a function of wavelength
INPUTS:
    T = Temperature of blackbody in K
OUTPUTS:
    Planck function as a function of wavelength (photons/m^2/s/m)
"""
#Planck function as a function of temperature and wavelength
def Planck(T,lam):
    return 2*np.pi*c/(lam**4)/(np.exp(h*c/(lam*k_B*T))-1)


def Planck_wrapper(T):
    def func(lam):
        return Planck(T,lam)
    return func


#Limb darkening as a function of the fraction of stellar radius (1 = stellar radius)
def limb_darkening(r):
    mu = np.sqrt(1-r**2)
    return 1-0.47*(1-mu)-0.23*(1-mu)**2

#How many multiples do I need of the baseline to keep the baseline between 10 and 600m
def baseline_checker(baseline):
    if baseline > 600:
        return 600, 0
    elif baseline >= 10:
        return baseline, 1
    else:
        n = 1
        while baseline < 10:
            n += 1
            baseline *= (2*n-1)
        return baseline, n

#Calculate zodiacal background in phot/s (no collecting area dependence)
def zodiacal_background(star,spec):

    #Get JWST background
    bg = jbt.background(star.RA, star.Dec, 10)

    #Get Zodi background
    year_zodi = bg.bkg_data['zodi_bg'][:,0]
    #index = np.argmin(year_zodi)
    index = np.argwhere(year_zodi == np.percentile(year_zodi, 30, interpolation='nearest'))[0,0]


    zodi_data = bg.bkg_data['zodi_bg'][index]*1e6*1e-26 #W/m^2/Hz/sr
    waves = bg.bkg_data['wave_array']*1e-6 #m

    #Interpolate the zodiacal background onto the spectrograph grid
    filter_waves = np.linspace(spec.wave_min,spec.wave_max,2000)

    f_zodi = interp1d(waves,zodi_data)

    #multiply by filter transmission
    zodi_radiance_hz = f_zodi(filter_waves) #W/m^2/Hz/sr

    #Convert to photons/s/m:
    # a) to wavelengths by multiplying by c/lam^2
    # b) to photons by dividing through by hc/lam
    # c) solid angle dependence: Omega \propto lam^2/Area
    # d) thus multiplying by Omega*Area = lam^2 removes solid angle and area dependence
    # Putting it together:
    zodi_spec_power = zodi_radiance_hz*filter_waves/h #phot/s/m

    #Sum over wavelength => Resultant power. Do it for the spectrograph channels
    zodi_power = []
    for power_split, wave_split in zip(np.split(zodi_spec_power,spec.n_channels),np.split(filter_waves,spec.n_channels)):
        zodi_power.append(np.trapz(power_split,wave_split)) #phot/s

    return np.array(zodi_power)


#Calculate planet signal flux (phot/s/m^2)
def calc_planet_signal(trans_map,planet,wave_pix2mas,spec,mode):

    signal_ls = []
    for flux,wave in zip(planet.flux,spec.channel_centres):
        planet_pos = planet.PAngSep/(wave_pix2mas*wave) #in pixels

        if mode == 1:
            p_trans = azimuthal_rms(trans_map,planet_pos)
        if mode == 2:
            p_trans = azimuthal_max(trans_map,planet_pos)

        p_flux = p_trans*flux
        signal_ls.append(p_flux)

    return np.array(signal_ls)


#Caclulate the local zodiacal radiance (photons/s/m^2/sr) along the zodiacal axis of symmetry
def calc_local_zodiacal_minimum(spec):

        #Ecliptic pole coordinates:
        ra = 270
        dec = 66.56

        #Get JWST background
        bg = jbt.background(ra, dec, 10)

        #Get Zodi background
        year_zodi = bg.bkg_data['zodi_bg'][:,0]
        index = np.argmin(year_zodi)

        zodi_data = bg.bkg_data['zodi_bg'][index]*1e6*1e-26 #W/m^2/Hz/sr
        waves = bg.bkg_data['wave_array']*1e-6 #m

        #Interpolate the zodiacal background onto the spectrograph grid
        filter_waves = np.linspace(spec.wave_min,spec.wave_max,2000)

        f_zodi = interp1d(waves,zodi_data)

        #multiply by filter transmission
        zodi_radiance_hz = f_zodi(filter_waves) #W/m^2/Hz/sr

        #Convert to photons/s/m^2/m/sr:
        # a) to wavelengths by multiplying by c/lam^2
        # b) to photons by dividing through by hc/lam
        # Putting it together:
        zodi_spec_rad = zodi_radiance_hz/(filter_waves*h) #photons/s/m^2/m/sr

        #Sum over wavelength => Resultant radiance. Do it for each spectral channel
        zodi_rad = []
        for rad_split, wave_split in zip(np.split(zodi_spec_rad,spec.n_channels),np.split(filter_waves,spec.n_channels)):
            #print(rad_split)
            #print(wave_split)
            zodi_rad.append(np.trapz(rad_split,wave_split)) #photons/s/m^2/sr

        return np.array(zodi_rad)


#Calc the exozodiacal flux (photons/s/m^2)
#Essentially sums over solid angle.
def calc_exozodiacal(star,trans_map,local_exozodi,pix2mas,sz,spec):

    #Create array
    arr = np.arange(sz)-sz/2
    x,y = np.meshgrid(arr,arr)
    r = np.sqrt(x**2 + y**2)

    r = r*pix2mas/rad2mas*star.Dist*pc/au #convert into au scale

    #Flux distribution from WISPR Observations (Stenborg 2020)
    # au = solar_rad/215
    r_in = 3/215 #3 Solar radii
    r_out = 19/215 #19 Solar radii
    lambda_r = np.piecewise(r, [r < r_in, ((r >= r_in) & (r <= r_out)), r > r_out],
                               [0, lambda x: (x-r_in)/(r_out-r_in), 1])

    #Column density distribution of zodiacal dust (IS THIS SCALING TRUE???)
    column_density = lambda_r*r**(-0.3)
    column_density[int(sz/2),int(sz/2)] = 0

    #Temperature distribution of zodiacal dust (IS THIS SCALING TRUE???)
    temp_dist = 300*r**(-0.5)
    temp_dist[int(sz/2),int(sz/2)] = 1

    #Solid angle of each pixel
    solid_angle_pix = (pix2mas/rad2mas)**2

    #Calculate exozodiacal flux for each wavelength channel
    exozodi = []
    for i in range(len(local_exozodi)):
        #Zodi flux at 1au is 2*local amount*exozodis
        local_scale_factor = 2*local_exozodi[i]*star.Exzod

        #Normalised planck distribution (i.e radiance)
        wavelength_sample = np.linspace(spec.channel_borders[i],spec.channel_borders[i]+spec.dlambda,50)

        planck_arr = []
        planck_norm = []
        for lam in wavelength_sample:
            planck_arr.append(Planck(temp_dist,lam))
            planck_norm.append(Planck(300,lam))

        planck_dist = np.trapz(np.array(planck_arr),wavelength_sample,axis=0)/np.trapz(np.array(planck_norm),wavelength_sample)

        flux_dist = column_density*planck_dist

        exozodi.append(np.sum(flux_dist*trans_map*local_scale_factor*solid_angle_pix))

    return np.array(exozodi)


def azimuthal_max(image,r):

    n_angles = 10000
    angles = np.linspace(0,2*np.pi,n_angles)

    centre = (int(image.shape[0]/2),int(image.shape[1]/2))

    #Planet out of field of view!
    if r > image.shape[0]/2:
        return 0

    max = 0
    for theta in angles:
        x = centre[0] + r*np.cos(theta)
        y = centre[1] + r*np.sin(theta)

        a = image[int(x),int(y)]

        if a > max:
            max = a

    return max

def azimuthal_rms(image,r):

    n_angles = 10000
    angles = np.linspace(0,2*np.pi,n_angles)

    centre = (int(image.shape[0]/2),int(image.shape[1]/2))

    #Planet out of field of view!
    if r > image.shape[0]/2:
        return 0

    sum = 0
    for theta in angles:
        x = centre[0] + r*np.cos(theta)
        y = centre[1] + r*np.sin(theta)

        a = image[int(x),int(y)]

        sum += a**2

    return np.sqrt(sum/n_angles)
