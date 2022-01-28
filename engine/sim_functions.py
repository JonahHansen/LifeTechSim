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


########################################################

"""
Spectrograph class that contains spectral parameters
Initiate units in um. Attributes in m

Inputs:
    wave_min = minimum wavelength
    wave_max = maximum wavelength
    baseline_wave = wave to optimise baselines around
    n_channels = number of spectral channels
"""
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
The spectral flux density for a blackbody
    T = Temperature of blackbody in K
    lam = wavelength in m
OUTPUTS:
    Spectral flux density
"""
def Planck(T,lam):
    with np.errstate(over='ignore', invalid='ignore'):
        return 2*np.pi*c/(lam**4)/(np.exp(h*c/(lam*k_B*T))-1)

"""
Return a function for the spectral flux density as a function of wavelength
Essentially a wrapper for the above function
INPUTS:
    T = Temperature of blackbody in K
OUTPUTS:
    Planck function as a function of wavelength (photons/m^2/s/m)
"""
def Planck_wrapper(T):
    def func(lam):
        return Planck(T,lam)
    return func


#Limb darkening as a function of the fraction of stellar radius (1 = stellar radius)
def limb_darkening(r):
    with np.errstate(divide='ignore', invalid='ignore'):
        mu = np.sqrt(1-r**2)
    return 1-0.05*(1-mu)-0.10*(1-mu)**2


##################### Azimuthal functions ###################################

"""
Calculate the array over azimuthal angles for a given response map and
radial position

Inputs:
    image = response map
    r = radial position to find the maximum over angles

Outputs:
    array of transmission over azimuthal angles
"""
def azimuthal_array(image,r):

    n_angles = 5000
    angles = np.linspace(0,2*np.pi,n_angles)

    centre = (int(image.shape[0]/2),int(image.shape[1]/2))

    #Planet out of field of view!
    if r > image.shape[0]/2:
        return 0

    arr = []
    for theta in angles:
        x = centre[0] + r*np.cos(theta)
        y = centre[1] + r*np.sin(theta)

        a = image[int(x),int(y)]

        arr.append(a)

    return np.array(arr)

"""
Calculate the RMS average over azimuthal angles for a given response map and
radial position

Inputs:
    image = map to average over
    r = radial position to average over angles

Outputs:
    RMS average over azimuth
"""
def azimuthal_rms(image,r):

    n_angles = 5000
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


"""
Calculate the mean average over azimuthal angles for a given response map and
radial position

Inputs:
    image = map to average over
    r = radial position to average over angles

Outputs:
    mean over azimuth
"""
def azimuthal_mean(image,r):

    n_angles = 5000
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

        sum += a

    return sum/n_angles

##################### Main Signal/Noise functions #########################

"""
Calculate planet signal flux (phot/s/m^2)

Inputs:
    outputs = list of (response,kernel) map tuples
    planet = planet (of planet class) to calculate the signal of
    wave_pix2mas = pixel to mas conversion factor, decoupled from wavelength
                   multiply by wavelength to get proper scale factor
    spec = spectrograph class that contains the spectral parameters
    mode = 1 = search, 2 = characterisation

Output:
    list of planet signal fluxes in phot/s/m^2 (for each kernel and wavelength)
"""
def calc_planet_signal(outputs,planet,wave_pix2mas,spec,mode):

    signal = []

    for (res,ker) in outputs:

        temp_signal = []
        for flux,wave in zip(planet.flux,spec.channel_centres):
            planet_pos = planet.PAngSep/(wave_pix2mas*wave) #in pixels

            if mode == 1:
                #Search mode is rms average of kernel azimuth
                p_trans = azimuthal_rms(ker,planet_pos)
                p_flux = p_trans*flux #multiply by flux
                temp_signal.append(p_flux)

            if mode == 2:
                #Characterisation mode is maximum of kernel azimuth
                p_flux_angle_array = np.abs(azimuthal_array(ker,planet_pos))*flux
                temp_signal.append(p_flux_angle_array)

        signal.append(np.array(temp_signal))

    if mode == 2:
        signal = np.array(signal)
        summed_signal = np.sum(signal,axis=(0,1))
        arg = np.argmax(summed_signal)
        signal = signal[:,:,arg]

    return signal


"""
Calculate planet shot flux (phot/s/m^2)

Inputs:
    outputs = list of (response,kernel) map tuples
    planet = planet (of planet class) to calculate the signal of
    wave_pix2mas = pixel to mas conversion factor, decoupled from wavelength
                   multiply by wavelength to get proper scale factor
    spec = spectrograph class that contains the spectral parameters
    mode = 1 = search, 2 = characterisation

Output:
    list of planet shot fluxes in phot/s/m^2 (for each kernel and wavelength)
"""
def calc_shot_noise(outputs,planet,wave_pix2mas,spec,mode):

    shot_noise = []

    for (res,ker) in outputs:

        temp_shot_noise = []
        for flux,wave in zip(planet.flux,spec.channel_centres):
            planet_pos = planet.PAngSep/(wave_pix2mas*wave) #in pixels

            if mode == 1:
                #Search mode is rms average of kernel azimuth
                p_trans = azimuthal_mean(res,planet_pos)
                p_flux = p_trans*flux #multiply by flux
                temp_shot_noise.append(p_flux)

            if mode == 2:
                #Characterisation mode is maximum of kernel azimuth
                p_flux_angle_array = np.abs(azimuthal_array(res,planet_pos))*flux
                temp_shot_noise.append(p_flux_angle_array)

        shot_noise.append(np.array(temp_shot_noise))


    if mode == 2:
        shot_noise = np.array(shot_noise)
        summed_shot = np.sum(shot_noise,axis=(0,1))
        arg = np.argmax(summed_shot)
        shot_noise = shot_noise[:,:,arg]

    return shot_noise


"""
Calculate stellar leakage flux (phot/s/m^2)
Uses an increased resolution simulation (like exozodiacal light)

Inputs:
    star = star (of star class) to calculate leakage of
    response_func = function (based on architecture) to calculate the list of response maps
    baseline = (shortest) baseline of the array
    base_wavelength = wavelength used to define the baseline

Output:
    list of stellar leakage fluxes in phot/s/m^2 (for each kernel and wavelength)
"""
def stellar_leakage(star,response_func,baseline,spec):

    sz = 500

    #calculate transmission (field of view has a radius equal to 5x star's angular radius)
    outputs = response_func(baseline,10*star.angRad/rad2mas,sz,spec.baseline_wave)

    #Create an array
    arr = np.arange(sz)-sz/2
    x,y = np.meshgrid(arr,arr)

    #Calculate leakage for each output
    leakage = []
    for (res,ker) in outputs:

        temp_leakage = []
        for i, wave in enumerate(spec.channel_centres):

            r = np.sqrt(x**2 + y**2)
            r *= 10/sz #Full grid is 5x the radius

            #Scale according to wavelength
            r *= wave/spec.baseline_wave

            #Limb_darkening, normalised over the total area
            I = limb_darkening(r)
            I[np.isnan(I)] = 0
            I/=np.sum(I) #normalise by total sum

            #Turn into limb_darkened flux

            I = I*star.flux[i]

            temp_leakage.append(np.sum(res*I))

        leakage.append(temp_leakage)

    return leakage


"""
Calculate zodiacal background power (phot/s) (no collecting area dependence)

Inputs:
    star = star (of star class) to calculate leakage of
    spec = spectrograph class that contains the spectral parameters

Output:
    list of zodiacal powers in phot/s (for each wavelength)
"""
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


"""
Calculate the minimum local zodiacal radiance (photons/s/m^2/sr) along the zodiacal axis of symmetry
That is, pointing in roughly the direction of the ecliptic pole

Inputs:
    spec = spectrograph class that contains the spectral parameters

Output:
    list of local zodiacal radiance in photons/s/m^2/sr (for each wavelength)
"""
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


"""
Calculate exozodiacal flux (phot/s/m^2)
Uses an resolution simulation, and assumes a distribution of the dust

Inputs:
    star = star (of star class) to calculate leakage of
    outputs = list of (response,kernel) map tuples
    local_zodi = list of minimum local zodiacal radiances (from above function)
    pix2mas = pixel to mas conversion for the default array map
    sz = size of response maps
    spec = spectrograph class that contains the spectral parameters

Output:
    list of exozodiacal fluxes in phot/s/m^2 (for each kernel and wavelength)
"""
#Calc the exozodiacal flux (photons/s/m^2)
#Essentially sums over solid angle.
def calc_exozodiacal(star,outputs,local_zodi,wave_pix2mas,sz,spec):

    #Create array
    arr = np.arange(sz)-sz/2
    x,y = np.meshgrid(arr,arr)

    #Calculate over all outputs
    exozodi = []
    for (res,ker) in outputs:

        #Calculate exozodiacal flux for each wavelength channel
        temp_exozodi = []
        for i, wave in enumerate(spec.channel_centres):

            pix2mas = wave_pix2mas*wave

            r = np.sqrt(x**2 + y**2)
            r = r*pix2mas/rad2mas*star.Dist*pc/au #convert into au scale

            #Flux distribution from WISPR Observations (Stenborg 2020)
            # au = solar_rad/215
            r_in = 3/215 #3 Solar radii
            r_out = 19/215 #19 Solar radii
            lambda_r = np.piecewise(r, [r < r_in, ((r >= r_in) & (r <= r_out)), r > r_out],
                                       [0, lambda x: (x-r_in)/(r_out-r_in), 1])

            #Column density distribution of zodiacal dust (IS THIS SCALING TRUE...) YES SEE BELOW
            #EXO-ZODI MODELLING FOR THE LARGE BINOCULAR TELESCOPE INTERFEROMETER - Kennedy 2015!
            with np.errstate(divide='ignore', invalid='ignore'):
                column_density = lambda_r*r**(-0.34)
            column_density[int(sz/2),int(sz/2)] = 0

            #Temperature distribution of zodiacal dust (IS THIS SCALING TRUE???)
            with np.errstate(divide='ignore', invalid='ignore'):
                temp_dist = 300*r**(-0.5)
            temp_dist[int(sz/2),int(sz/2)] = 1

            #Solid angle of each pixel
            solid_angle_pix = (pix2mas/rad2mas)**2

            #Zodi flux at 1au is 2*local amount*exozodis
            local_scale_factor = 2*local_zodi[i]*star.Exzod

            #Sample each wavelength channeel
            wavelength_sample = np.linspace(spec.channel_borders[i],spec.channel_borders[i]+spec.dlambda,5)

            planck_arr = []
            planck_norm = []
            #Calculate planck functions for each wavelength samples
            for lam in wavelength_sample:
                #Calculate the planck function for the temperature distribution
                planck_arr.append(Planck(temp_dist,lam))
                #Calculate the planck function for 300K (to normalise by)
                planck_norm.append(Planck(300,lam))

            #Normalised planck distribution (i.e radiance)
            planck_dist = np.trapz(np.array(planck_arr),wavelength_sample,axis=0)/np.trapz(np.array(planck_norm),wavelength_sample)

            #Multiply by column density for flux distribution
            flux_dist = column_density*planck_dist

            #Exozodiacal is the sum over the array of the flux distribution, transmission,
            #local scale factor and the pixel solid angle
            temp_exozodi.append(np.sum(flux_dist*res*local_scale_factor*solid_angle_pix))
        exozodi.append(temp_exozodi)

    return exozodi
