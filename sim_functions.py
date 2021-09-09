import numpy as np
from jwst_backgrounds import jbt
from scipy.interpolate import interp1d


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


def zodiacal_background(star,filter):

    #Get JWST background
    bg = jbt.background(star.RA, star.Dec, 10)

    #Get Zodi background
    min_index = np.argmin(bg.bkg_data['zodi_bg'][:,0])
    zodi_data = bg.bkg_data['zodi_bg'][min_index]*1e6*1e-26 #W/m^2/Hz/sr
    waves = bg.bkg_data['wave_array']*1e-6 #m

    #Interpolate the filter and the zodiacal background
    f_filter = interp1d(filter.Wavel,filter.Trans)
    f_zodi = interp1d(waves,zodi_data)

    #Common wavelengths
    common_waves = np.linspace(filter.Wavel[0],filter.Wavel[1],1000)

    #Integrate over solid angle and multiply by filter transmission
    zodi_irradiance_hz = np.pi*f_filter(common_waves)*f_zodi(common_waves) #W/m^2/Hz
    #Convert to photons ((c/lam**2)/(hc/lam))
    zodi_irradiance_phot = zodi_irradiance_hz/(h*common_waves) #phot/m^2/s/m

    #Sum over wavelength
    zodi_flux = np.trapz(zodi_irradiance_m,common_waves) #phot/m^2/s

    return zodi_flux


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
