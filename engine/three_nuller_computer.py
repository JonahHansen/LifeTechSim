import sys
sys.path.append("..")
import numpy as np
from astropy import constants as const
import engine.sim_functions as sf
from opticstools import knull
from opticstools.opticstools import azimuthalAverage

rad2mas = np.degrees(1)*3600e3 #Number of milliarcsec in one radian

#The baseline given to this function is the one that defines the "optimised" transmission map
#The baseline is the distance between adjacent spacecraft
def triangle(baseline):
    angles = np.linspace(0,2*np.pi,4)
    xs = 0.5773*baseline*np.sin(angles)
    ys = 0.5773*baseline*np.cos(angles)
    return np.array([xs,ys]).T[:-1]


def get_nuller_response(baseline,fov,sz,base_wavelength):

    N = knull.make_nuller_mat3()

    telescope_array = triangle(baseline)

    sky_angles = np.linspace(-fov/2,fov/2,sz)

    xy = np.meshgrid(sky_angles, sky_angles, indexing='ij')

    #x,y are telescope positions in units of wavelength
    x = telescope_array[:,0]/base_wavelength
    y = telescope_array[:,1]/base_wavelength

    #Response is the 5 output electric fields as a function of the position on the sky
    response = np.zeros((3,sz,sz), dtype='complex')

    for i in range(3):
        for k in range(3):
            #Inputs have a phase equal to xy array - linear multiplied by spatial frequency
            #ul + vm?
            response[k] += np.exp(2*np.pi*1j*(xy[0]*x[i] + xy[1]*y[i]))*N[k,i] #

    response = np.abs(response)**2
    response /= (np.max(response[0])/3) #normalise by flux per telescope

    #Create the kernel nulls. This is turning the output intensities into the kernel nulls (K in 2018 paper)
    k = response[1]-response[2]

    return response, k #return intensity per telescope


#Calculate stellar leakage flux through an increased resolution simulation (like exozodiacal light)
def stellar_leakage(star,baseline,base_wavelength):

    sz = 400

    #Create an array
    arr = np.arange(sz)-sz/2
    x,y = np.meshgrid(arr,arr)
    r = np.sqrt(x**2 + y**2)
    r *= 2/sz #convert into multiples of stellar radius

    pixel_size = r[100,100] - r[100,99]

    #Limb_darkening, normalised over the total area
    I = sf.limb_darkening(r)#
    I[np.isnan(I)] = 0
    I/=np.sum(I)

    #Turn into limb_darkened flux
    fluxes = np.tile(star.flux,(sz,sz,1)).T
    I = I*fluxes

    #calculate transmission (field of view has a radius equal to star's angular radius)
    res, k = get_nuller_response(baseline,2*star.angRad/rad2mas,sz,base_wavelength)

    #Calculate leakage
    leakage = np.sum(res[1]*I,axis=(1,2))

    return leakage


#Calculate stellar leakage flux through Mike's method
def Mike_stellar_leakage(star,response,pix2mas):
    #Find the null function numerically
    #Integral over theta, divide by two pi, as a function of radius
    r_ix,y2 = azimuthalAverage(response[1], returnradii=True, binsize=0.8)
    #averaging response function over theta as a function of radius. yi is average value

    #Dodgy fit for a quartic and parabola
    second_order_coeff = np.median(y2[1:16]/r_ix[1:16]**2)/pix2mas**2

    r = np.linspace(0,1,200)
    I = sf.limb_darkening(r)

    #Total leakage is an integral of coeff multiplied by r^2 or r^4 and limb darkening intensity
    mn_r2 = (np.trapz(I*r**3, r)/np.trapz(I*r, r))**.5

    #convert to flux...
    leakage = second_order_coeff*(mn_r2*star.angRad)**2*star.flux
    return leakage


def compute(star,mode,spec,sz,scale_factor,local_exozodi):

    #Define wavelength to fix baseline
    base_wavelength = spec.baseline_wave #

    print("\nCalculating Zodiacal")

    #zodiacal power (phot/s) per telescope
    zodiacal = sf.zodiacal_background(star,spec)

    print("\nCalculating Response")

    ls_row_data = []

    if mode == 1:
        baseline = base_wavelength/2*rad2mas/star.HZAngle
        fov = base_wavelength/baseline*scale_factor
        response, k = get_nuller_response(baseline,fov,sz,base_wavelength)
        pix2mas = fov*rad2mas/sz

        print("\nCalculating Exozodiacal")
        #exozodiacal flux (phot/s/m^2) per telescope
        exozodiacal = sf.calc_exozodiacal(star,response[1],local_exozodi,pix2mas,sz,spec)

        print("\nCalculating Leakage")
        #Calc stellar leakage flux (phot/s/m^2) per telescope
        leakage = stellar_leakage(star,baseline,base_wavelength)

        """
        #Calc stellar leakage flux (phot/s/m^2) per telescope
        leakage = Mike_stellar_leakage(star,response,pix2mas)
        """

        for planet in star.Planets:

            #pix2mas conversion, removing wavelength dependence
            #multiply by wavelength to get conversion factor for that wavelength
            wave_pix2mas = pix2mas/base_wavelength

            print("\nCalculating Signal")
            #signal flux (phot/s/m^2) per telescope
            signal = sf.calc_planet_signal(k,planet,wave_pix2mas,spec,mode)

            row_data = {"star_name":star.Name, "planet_name":planet.Name, "universe_no":planet.UNumber,
                        "star_no":star.SNumber,"planet_no":planet.PNumber,"baseline (m)":baseline,
                        "array_angle (mas)":star.HZAngle, "planet_angle (mas)":planet.PAngSep, "star_type":star.Stype,
                        "star_flux (ph/s/m2)":star.flux,"planet_flux (ph/s/m2)":planet.flux,
                        "planet_temp (K)":planet.PTemp,"planet_radius (Earth_Rad)":planet.PRad,
                        "signal (ph/s/m2)":signal,"leakage (ph/s/m2)":leakage,"exozodiacal (ph/s/m2)":exozodiacal,
                        "zodiacal (ph/s)":zodiacal}

            ls_row_data.append(row_data)


    if mode == 2:
        for planet in star.Planets:
            baseline = base_wavelength/2*rad2mas/planet.PAngSep
            fov = base_wavelength/baseline*scale_factor
            response, k = get_nuller_response(baseline,fov,sz,base_wavelength)
            pix2mas = fov*rad2mas/sz

            print("\nCalculating Exozodiacal")
            #exozodiacal flux (phot/s/m^2) per telescope
            exozodiacal = sf.calc_exozodiacal(star,response[1],local_exozodi,pix2mas,sz,spec)

            print("\nCalculating Leakage")
            #Calc stellar leakage flux (phot/s/m^2) per telescope
            leakage = stellar_leakage(star,baseline,base_wavelength,sz)

            """
            #Calc stellar leakage flux (phot/s/m^2) per telescope
            leakage = Mike_stellar_leakage(star,response,pix2mas)
            """
            #pix2mas conversion, removing wavelength dependence
            #multiply by wavelength to get conversion factor for that wavelength
            wave_pix2mas = pix2mas/base_wavelength

            print("\nCalculating Signal")
            #signal flux (phot/s/m^2) per telescope
            signal = sf.calc_planet_signal(k,planet,wave_pix2mas,spec,mode)

            row_data = {"star_name":star.Name, "planet_name":planet.Name, "universe_no":planet.UNumber,
                        "star_no":star.SNumber,"planet_no":planet.PNumber,"baseline (m)":baseline,
                        "array_angle (mas)":planet.PAngSep, "planet_angle (mas)":planet.PAngSep, "star_type":star.Stype,
                        "star_flux (ph/s/m2)":star.flux,"planet_flux (ph/s/m2)":planet.flux,
                        "planet_temp (K)":planet.PTemp,"planet_radius (Earth_Rad)":planet.PRad,
                        "signal (ph/s/m2)":signal,"leakage (ph/s/m2)":leakage,"exozodiacal (ph/s/m2)":exozodiacal,
                        "zodiacal (ph/s)":zodiacal}

            ls_row_data.append(row_data)

    return ls_row_data
