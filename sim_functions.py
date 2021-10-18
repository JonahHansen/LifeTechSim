import numpy as np
from jwst_backgrounds import jbt
from scipy.interpolate import interp1d
from astropy import constants as const

h = const.h.value #Planck constant
au = const.au.value #AU in m
pc = const.pc.value #parsec in m
rad2mas = np.degrees(1)*3600e3 #Number of milliarcsec in one radian

#Initiate units in um. Attributes in m
class Spectrograph():
    def __init__(self, mean_wave, bandwidth, baseline_wave, n_channels):

        self.bandwidth = bandwidth*1e-6
        self.mean = mean_wave*1e-6
        self.wave_min = self.mean - self.bandwidth/2
        self.wave_max = self.mean + self.bandwidth/2

        self.baseline_wave = baseline_wave*1e-6

        self.channel_borders = np.linspace(self.wave_min,self.wave_max,n_channels+1)[:-1]
        self.dlambda = self.channel_borders[1]-self.channel_borders[0]
        self.channel_centres = (self.channel_borders + self.dlambda/2)

        self.eff_resolution = self.mean/self.dlambda



def limb_darkening(r):
    mu = np.sqrt(1-r**2)
    return 1-0.47*(1-mu)-0.23*(1-mu)**2


def get_SNR(signal,leakage,zodiacal,N,area,exp_time,eta):

    C_per_t = area*exp_time*eta
    #Signal
    signal_per_t = np.sum(signal*C_per_t)
    #Noise contributions: shot, leakage and zodiacal background
    noise_per_t = np.sqrt(np.sum(signal*C_per_t) + np.sum(leakage*C_per_t) + zodiacal*C_per_t)
    SNR_per_t = signal_per_t/noise_per_t


    C_tot = N*area*exp_time*eta
    #Signal
    signal_tot = np.sum(signal*C_tot)
    #Noise contributions: shot, leakage and zodiacal background
    noise_tot = np.sqrt(np.sum(signal*C_tot) + np.sum(leakage*C_tot) + zodiacal*C_tot)
    SNR_tot = signal_tot/noise_tot

    return SNR_per_t, SNR_tot


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

    #Sum over wavelength => Resultant power
    zodi_power = np.trapz(zodi_spec_power,filter_waves) #phot/s

    return zodi_power


#Calculate planet signal flux (phot/s/m^2)
def calc_planet_signal(trans_map,planet,wave_pix2mas):

    total = 0
    for flux,wave in zip(planet.flux,planet.flux_wave):
        planet_pos = planet.PAngSep/(wave_pix2mas*wave) #in pixels
        p_trans = azimuthal_rms(trans_map,planet_pos)
        p_flux = p_trans*flux
        total += p_flux

    return total


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

        #Sum over wavelength => Resultant power
        zodi_rad = np.trapz(zodi_spec_rad,filter_waves) #photons/s/m^2/sr

        return zodi_rad


#Calc the exozodiacal flux (photons/s/m^2)
#Essentially sums over solid angle.
def calc_exozodiacal(star,trans_map,local_exozodi,pix2mas,sz):

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
    flux_dist = lambda_r*r**(-2.3)
    flux_dist[int(sz/2),int(sz/2)] = 0

    #Zodi flux at 1au is 2*local amount*exozodis
    local_scale_factor = 2*local_exozodi*star.Exzod

    solid_angle_pix = (pix2mas/rad2mas)**2

    return np.sum(flux_dist*trans_map*local_scale_factor*solid_angle_pix)


def azimuthal_rms(image,r):
    n_angles = 10000
    angles = np.linspace(0,2*np.pi,n_angles)

    centre = (int(image.shape[0]/2),int(image.shape[1]/2))
    sum = 0
    for theta in angles:
        x = centre[0] + r*np.cos(theta)
        y = centre[1] + r*np.sin(theta)

        a = image[int(x),int(y)]

        sum += a**2

    return np.sqrt(sum/n_angles)
