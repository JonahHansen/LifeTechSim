import sys
sys.path.append("..")
import numpy as np
import engine.sim_functions as sf

rad2mas = np.degrees(1)*3600e3 #Number of milliarcsec in one radian

"""
Compute the relevant signal/noise fluxes for each planet around a given star and store in a dictionary

Inputs:
    star = star (of star class) to compute fluxes for
    mode = what mode of the interferometer? 1 = search, 2 = characterisation
    nuller_response = function to calculate the response maps. Changes depending on what architecture is used
    base_scale_factor = factor to scale the baseline. Changed based on architecture to optimise the transmission for a given sky angle
    fov_scale_factor = how large is the field of view? Given as a multiple of the optimised sky angle for a half FOV
    local_exozodi = the radiance of the local zodiacal background, to calculate the exozodiacal light

Output:
    dictionary of relevant fluxes and data, to be used in further analysis
"""
def compute(star,mode,nuller_response,spec,sz,base_scale_factor,fov_scale_factor,local_exozodi):

    #Define wavelength to fix baseline
    base_wavelength = spec.baseline_wave #

    #zodiacal background power (phot/s)
    zodiacal = sf.zodiacal_background(star,spec)

    ls_row_data = []

    #SEARCH MODE: array stays in same position per star
    if mode == 1:
        #Baseline is scaled by an appropriate factor for each architecture to maximise
        #transmission in the habitable zone
        baseline = base_scale_factor*base_wavelength*rad2mas/star.HZAngle

        #Half field of view is set as the scale_factor times the HZ angle
        fov = 2*fov_scale_factor*star.HZAngle/rad2mas

        #Get response maps
        outputs = nuller_response(baseline,fov,sz,base_wavelength)
        pix2mas = fov*rad2mas/sz #number of mas per pixel

        #pix2mas conversion, removing wavelength dependence
        #multiply by wavelength to get conversion factor for that wavelength
        wave_pix2mas = pix2mas/base_wavelength

        #exozodiacal flux (phot/s/m^2) per kernel
        exozodiacal = sf.calc_exozodiacal(star,outputs,local_exozodi,wave_pix2mas,sz,spec)

        #Calc stellar leakage flux (phot/s/m^2) per kernel
        leakage = sf.stellar_leakage(star,nuller_response,baseline,spec)

        for planet in star.Planets:

            print("\nCalculating Signal")
            #signal flux (phot/s/m^2) per kernel
            signal = sf.calc_planet_signal(outputs,planet,wave_pix2mas,spec,mode)

            #shot noise (phot/s/m^2) per kernel
            shot_noise = sf.calc_shot_noise(outputs,planet,wave_pix2mas,spec,mode)

            row_data = {"star_name":star.Name, "planet_name":planet.Name,
                        "universe_no":planet.UNumber,"star_no":star.SNumber,"planet_no":planet.PNumber,
                        "star_type":star.Stype,"star_distance":star.Dist,"baseline":baseline,
                        "array_angle":star.HZAngle, "planet_angle":planet.PAngSep,
                        "star_flux":star.flux,"planet_flux":planet.flux,
                        "planet_temp":planet.PTemp,"planet_radius":planet.PRad,
                        "habitable":str(planet.isHZ),
                        "signal":signal,
                        "shot":shot_noise,
                        "leakage":leakage,
                        "exozodiacal":exozodiacal,
                        "zodiacal":zodiacal}

            ls_row_data.append(row_data)

    #CHARACTERISATION MODE: array changes position for each planet
    if mode == 2:
        for planet in star.Planets:
            #Baseline is scaled by an appropriate factor for each architecture to maximise
            #transmission for the planet
            baseline = base_scale_factor*base_wavelength*rad2mas/planet.PAngSep

            #Half field of view is set as the scale_factor times the planet angular separation
            fov = 2*fov_scale_factor*planet.PAngSep/rad2mas

            #Get response maps
            outputs = nuller_response(baseline,fov,sz,base_wavelength)
            pix2mas = fov*rad2mas/sz #number of mas per pixel

            #pix2mas conversion, removing wavelength dependence
            #multiply by wavelength to get conversion factor for that wavelength
            wave_pix2mas = pix2mas/base_wavelength

            #exozodiacal flux (phot/s/m^2) per telescope
            exozodiacal = sf.calc_exozodiacal(star,outputs,local_exozodi,wave_pix2mas,sz,spec)

            #Calc stellar leakage flux (phot/s/m^2) per telescope
            leakage = sf.stellar_leakage(star,nuller_response,baseline,spec)

            #signal flux (phot/s/m^2) per telescope
            signal = sf.calc_planet_signal(outputs,planet,wave_pix2mas,spec,mode)

            #shot noise (phot/s/m^2) per telescope
            shot_noise = sf.calc_shot_noise(outputs,planet,wave_pix2mas,spec,mode)

            row_data = {"star_name":star.Name, "planet_name":planet.Name,
                        "universe_no":planet.UNumber,"star_no":star.SNumber,"planet_no":planet.PNumber,
                        "star_type":star.Stype,"star_distance":star.Dist,"baseline":baseline,
                        "array_angle":planet.PAngSep, "planet_angle":planet.PAngSep,
                        "star_flux":star.flux,"planet_flux":planet.flux,
                        "planet_temp":planet.PTemp,"planet_radius":planet.PRad,
                        "habitable":str(planet.isHZ),
                        "signal":signal,
                        "shot":shot_noise,
                        "leakage":leakage,
                        "exozodiacal":exozodiacal,
                        "zodiacal":zodiacal}

            ls_row_data.append(row_data)

    return ls_row_data
