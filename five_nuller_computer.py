import numpy as np
from astropy import constants as const
from sim_functions import limb_darkening, zodiacal_background, calc_exozodiacal, calc_planet_signal
from opticstools import knull
from opticstools.opticstools import azimuthalAverage

rad2mas = np.degrees(1)*3600e3 #Number of milliarcsec in one radian


def pentagon(baseline):
    angles = np.linspace(0,2*np.pi,6)
    xs = 1.176*baseline*np.sin(angles)
    ys = 1.176*baseline*np.cos(angles)
    return np.array([xs,ys]).T[:-1]


def get_nuller_response(baseline,fov,sz,base_wavelength):

    N = knull.make_nuller_mat5()

    telescope_array = pentagon(baseline)

    sky_angles = np.linspace(-fov/2,fov/2,sz)

    xy = np.meshgrid(sky_angles, sky_angles, indexing='ij')

    #x,y are telescope positions in units of wavelength
    x = telescope_array[:,0]/base_wavelength
    y = telescope_array[:,1]/base_wavelength

    #Response is the 5 output electric fields as a function of the position on the sky
    response = np.zeros((5,sz,sz), dtype='complex')

    for i in range(5):
        for k in range(5):
            #Inputs have a phase equal to xy array - linear multiplied by spatial frequency
            #ul + vm?
            response[k] += np.exp(2*np.pi*1j*(xy[0]*x[i] + xy[1]*y[i]))*N[k,i] #

    response = np.abs(response)**2
    response /= (np.max(response[0])/5) #normalise by flux per telescope

    #Create the kernel nulls. This is turning the output intensities into the kernel nulls (K in 2018 paper)

    k1 = response[2]-response[3]
    k2 = response[1]-response[4]

    return response, k1, k2 #return intensity per telescope


#Calculate stellar leakage flux through an increased resolution simulation (like exozodiacal light)
def stellar_leakage(star,baseline,base_wavelength,sz):

    arr = np.arange(sz)-sz/2
    x,y = np.meshgrid(arr,arr)
    r = np.sqrt(x**2 + y**2)
    r *= 2/sz #convert into multiples of stellar radius

    pixel_size = r[100,100] - r[100,99]

    #Scale the flux per pixel (Radius of 1 => star area of pi)
    Flux_per_pixel = star.Flux/np.pi*pixel_size**2

    I = limb_darkening(r)*Flux_per_pixel#limb darkened flux in phot/s/m^2

    I[np.isnan(I)] = 0

    #calculate transmission (field of view has a radius equal to star's angular radius)
    res, k1, k2 = get_nuller_response(baseline,2*star.angRad/rad2mas,sz,base_wavelength)

    #import pdb; pdb.set_trace()

    leakage_2nd = np.sum(res[1]*I)
    leakage_4th = np.sum(res[2]*I)

    return leakage_2nd, leakage_4th


#Calculate stellar leakage flux through Mike's method
def Mike_stellar_leakage(star,response,pix2mas):
    #Find the null function numerically
    #Integral over theta, divide by two pi, as a function of radius
    r_ix,y2 = azimuthalAverage(response[1], returnradii=True, binsize=0.8)
    r_ix,y4 = azimuthalAverage(response[2], returnradii=True, binsize=0.8)
    #averaging response function over theta as a function of radius. yi is average value

    #Dodgy fit for a quartic and parabola
    second_order_coeff = np.median(y2[1:16]/r_ix[1:16]**2)/pix2mas**2
    fourth_order_coeff = np.median(y4[1:16]/r_ix[1:16]**4)/pix2mas**4

    r = np.linspace(0,1,200)
    I = limb_darkening(r)

    #Total leakage is an integral of coeff multiplied by r^2 or r^4 and limb darkening intensity
    mn_r2 = (np.trapz(I*r**3, r)/np.trapz(I*r, r))**.5
    mn_r4 = (np.trapz(I*r**5, r)/np.trapz(I*r, r))**.25

    #convert to flux...
    leakage_2nd = second_order_coeff*(mn_r2*star.angRad)**2*star.Flux
    leakage_4th = fourth_order_coeff*(mn_r4*star.angRad)**4*star.Flux

    return leakage_2nd, leakage_4th


def compute(star,mode,spec,sz,scale_factor,local_exozodi):

    #Define wavelength to fix baseline
    base_wavelength = spec.baseline_wave #

    print("\nCalculating Zodiacal")

    #zodiacal power (phot/s) per telescope
    zodiacal = zodiacal_background(star,spec)

    print("\nCalculating Response")

    ls_row_data = []

    if mode == 1:
        baseline = base_wavelength/2*rad2mas/star.HZAngle
        fov = base_wavelength/baseline*scale_factor
        response, k1, k2 = get_nuller_response(baseline,fov,sz,base_wavelength)
        pix2mas = fov*rad2mas/sz


        print("\nCalculating Exozodiacal")
        #exozodiacal flux (phot/s/m^2) per telescope
        ave_trans_map = 0.5*(response[1] + response[2])
        exozodiacal = calc_exozodiacal(star,ave_trans_map,local_exozodi,pix2mas,sz)


        print("\nCalculating Leakage")

        #Calc stellar leakage flux (phot/s/m^2) per telescope
        leakage_2nd,leakage_4th = stellar_leakage(star,baseline,base_wavelength,sz)

        """
        #Calc stellar leakage flux (phot/s/m^2) per telescope
        leakage_2nd,leakage_4th = Mike_stellar_leakage(star,response,pix2mas)
        """

        for planet in star.Planets:

            #pix2mas conversion, removing wavelength dependence
            #multiply by wavelength to get conversion factor for that wavelength
            wave_pix2mas = pix2mas/base_wavelength

            print("\nCalculating Signal")
            #signal flux (phot/s/m^2) per telescope
            signal_k1 = calc_planet_signal(k1,planet,wave_pix2mas)
            signal_k2 = calc_planet_signal(k2,planet,wave_pix2mas)

            row_data = {"star_name":star.Name, "planet_name":planet.Name, "universe_no":planet.UNumber,
                        "star_no":star.SNumber,"planet_no":planet.PNumber,
                        "signal_k1":signal_k1,"leakage_2nd":leakage_2nd,
                        "signal_k2":signal_k2,"leakage_4th":leakage_4th,
                        "zodiacal":zodiacal,"exozodiacal":exozodiacal}

            ls_row_data.append(row_data)


    if mode == 2:
        for planet in star.Planets:
            baseline = base_wavelength/2*rad2mas/planet.PAngSep
            fov = base_wavelength/baseline*scale_factor
            response, k1, k2 = get_nuller_response(baseline,fov,sz,base_wavelength)
            pix2mas = fov*rad2mas/sz

            print("\nCalculating Exozodiacal")

            ave_trans_map = 0.5*(response[1] + response[2])

            #exozodiacal flux (phot/s/m^2) per telescope
            exozodiacal = calc_exozodiacal(star,ave_trans_map,local_exozodi,pix2mas,sz)

            print("\nCalculating Leakage")
            """
            #Calc stellar leakage flux (phot/s/m^2) per telescope
            leakage_2nd,leakage_4th = stellar_leakage(star,baseline,base_wavelength,sz)

            """
            #Calc stellar leakage flux (phot/s/m^2) per telescope
            leakage_2nd,leakage_4th = Mike_stellar_leakage(star,response,pix2mas)


            #pix2mas conversion, removing wavelength dependence
            #multiply by wavelength to get conversion factor for that wavelength
            wave_pix2mas = pix2mas/base_wavelength

            print("\nCalculating Signal")
            #signal flux (phot/s/m^2) per telescope
            signal_k1 = calc_planet_signal(k1,planet,wave_pix2mas)
            signal_k2 = calc_planet_signal(k2,planet,wave_pix2mas)

            row_data = {"star_name":star.Name, "planet_name":planet.Name, "universe_no":planet.UNumber,
                        "star_no":star.SNumber,"planet_no":planet.PNumber,
                        "signal_k1":signal_k1,"leakage_2nd":leakage_2nd,
                        "signal_k2":signal_k2,"leakage_4th":leakage_4th,
                        "zodiacal":zodiacal,"exozodiacal":exozodiacal}

            ls_row_data.append(row_data)

    return ls_row_data
