import numpy as np
from astropy import constants as const
from sim_functions import limb_darkening, get_SNR, zodiacal_background
from opticstools import knull

rad2mas = np.degrees(1)*3600e3 #Number of milliarcsec in one radian


def pentagon(baseline):
    angles = np.linspace(0,2*np.pi,6)
    xs = 1.176*baseline*np.sin(angles)
    ys = 1.176*baseline*np.cos(angles)
    return np.array([xs,ys]).T[:-1]


def compute(star,mode,filter,sz,scale_factor,local_exozodi):

    #Define wavelength to fix baseline
    wavelength = filter.mean #

    #zodiacal power (phot/s) per telescope
    zodiacal = zodiacal_background(star,filter)

    def get_nuller_response(baseline):

        N = knull.make_nuller_mat5()

        telescope_array = pentagon(baseline)

        fov = wavelength/baseline*scale_factor
        sky_angles = np.linspace(-fov/2,fov/2,sz)

        xy = np.meshgrid(sky_angles, sky_angles, indexing='ij')

        #x,y are telescope positions in units of wavelength
        x = telescope_array[:,0]/wavelength
        y = telescope_array[:,1]/wavelength

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

    ls_row_data = []

    if mode == 1:
        baseline = wavelength/2*rad2mas/star.HZAngle
        response, k1, k2 = get_nuller_response(baseline)
        pix2mas = wavelength/baseline*rad2mas*scale_factor/sz

        #exozodiacal flux (phot/s/m^2) per telescope
        exozodiacal = calc_exozodiacal(star,response[1],response[2],local_exozodi,pix2mas,sz)

        #Calc stellar leakage transmission
        leakage = s_trans*star.Flux

        for planet in star.Planets:

            #signal flux (phot/s/m^2) per telescope
            signal_k1 = calc_planet_signal(k1,planet,pix2mas)
            signal_k2 = calc_planet_signal(k2,planet,pix2mas)

            row_data = {"star_name":star.SName, "planet_name":planet.Name, "universe_no":planet.UNumber,
                        "star_no":star.SNumber,"planet_no":planet.PNumber,
                        "signal":signal,"leakage":leakage,"zodiacal":zodiacal,
                        "exozodiacal":exozodiacal}

            ls_row_data.append(row_data)


    if mode == 2:
        for planet in star.Planets:
            baseline = wavelength/2*rad2mas/planet.PAngSep
            response, k1, k2 = get_nuller_response(baseline)
            pix2mas = wavelength/baseline*rad2mas*scale_factor/sz

            #exozodiacal flux (phot/s/m^2) per telescope
            exozodiacal = calc_exozodiacal(star,response[1],response[2],local_exozodi,pix2mas,sz)

            #Calc stellar leakage transmission
            leakage = s_trans*star.Flux

            #signal flux (phot/s/m^2) per telescope
            signal_k1 = calc_planet_signal(k1,planet,pix2mas)
            signal_k2 = calc_planet_signal(k2,planet,pix2mas)

            row_data = {"star_name":star.SName, "planet_name":planet.Name, "universe_no":planet.UNumber,
                        "star_no":star.SNumber,"planet_no":planet.PNumber,
                        "signal_k1":signal_k1,"leakage_k1":leakage_k1,
                        "signal_k2":signal_k2,"leakage_k2":leakage_k2,
                        "zodiacal":zodiacal,"exozodiacal":exozodiacal}

            ls_row_data.append(row_data)

    return ls_row_data


#Find the null function numerically
#Integral over theta, divide by two pi, as a function of radius
r_ix,y2 = azimuthalAverage(response[1], returnradii=True, binsize=0.8)
r_ix,y4 = azimuthalAverage(response[2], returnradii=True, binsize=0.8)
#averaging response function over theta as a function of radius. yi is average value


#Dodgy fit for a quartic and parabola
second_order_coeff = np.median(y2[1:16]/r_ix[1:16]**2)/mas_pix**2
fourth_order_coeff = np.median(y4[1:16]/r_ix[1:16]**2)/mas_pix**2


r = np.linspace(0,1,200)

#Total leakage is an integral of coeff multiplied by r^2 or r^4 and limb darkening intensity
mn_r2 = (np.trapz(I*r**3, r)/np.trapz(I*r, r))**.5
mn_r4 = (np.trapz(I*r**5, r)/np.trapz(I*r, r))**.25

#tau Ceti calculation (Scaling by radius of star)
star_mas = Rstar/214*plx
leakage_2nd = second_order_coeff*(mn_r2*star_mas)**2
leakage_4th = fourth_order_coeff*(mn_r4*star_mas)**4
