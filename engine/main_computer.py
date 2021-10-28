import sys
sys.path.append("..")
import numpy as np
import engine.sim_functions as sf

rad2mas = np.degrees(1)*3600e3 #Number of milliarcsec in one radian


def compute(star,mode,nuller_response,spec,sz,base_scale_factor,fov_scale_factor,local_exozodi):

    #Define wavelength to fix baseline
    base_wavelength = spec.baseline_wave #

    print("\nCalculating Zodiacal")

    #zodiacal power (phot/s) per telescope
    zodiacal = sf.zodiacal_background(star,spec)

    print("\nCalculating Response")

    ls_row_data = []

    if mode == 1:
        baseline = base_scale_factor*base_wavelength*rad2mas/star.HZAngle

        #INSERT BASELINE CHECKER BASED ON OPTIMISATION?

        fov = 2*fov_scale_factor*star.HZAngle/rad2mas
        outputs = nuller_response(baseline,fov,sz,base_wavelength)
        pix2mas = fov*rad2mas/sz

        print("\nCalculating Exozodiacal")
        #exozodiacal flux (phot/s/m^2) per telescope
        exozodiacal = sf.calc_exozodiacal(star,outputs,local_exozodi,pix2mas,sz,spec)

        print("\nCalculating Leakage")
        #Calc stellar leakage flux (phot/s/m^2) per telescope
        leakage = sf.stellar_leakage(star,nuller_response,baseline,base_wavelength)

        for planet in star.Planets:

            #pix2mas conversion, removing wavelength dependence
            #multiply by wavelength to get conversion factor for that wavelength
            wave_pix2mas = pix2mas/base_wavelength

            print("\nCalculating Signal")
            #signal flux (phot/s/m^2) per telescope
            signal = sf.calc_planet_signal(outputs,planet,wave_pix2mas,spec,mode)

            shot_noise = sf.calc_shot_noise(outputs,planet,wave_pix2mas,spec,mode)

            row_data = {"star_name":star.Name, "planet_name":planet.Name,
                        "universe_no":planet.UNumber,"star_no":star.SNumber,"planet_no":planet.PNumber,
                        "star_type":star.Stype,"star_distance (pc)":star.Dist,"baseline (m)":baseline,
                        "array_angle (mas)":star.HZAngle, "planet_angle (mas)":planet.PAngSep,
                        "star_flux (ph/s/m2)":star.flux,"planet_flux (ph/s/m2)":planet.flux,
                        "planet_temp (K)":planet.PTemp,"planet_radius (Earth_Rad)":planet.PRad,
                        "signal (ph/s/m2)":signal,
                        "shot (ph/s/m2)":shot_noise,
                        "leakage (ph/s/m2)":leakage,
                        "exozodiacal (ph/s/m2)":exozodiacal,
                        "zodiacal (ph/s)":zodiacal}

            ls_row_data.append(row_data)


    if mode == 2:
        for planet in star.Planets:
            baseline = base_scale_factor*base_wavelength*rad2mas/planet.PAngSep
            fov = 2*fov_scale_factor*planet.PAngSep/rad2mas
            outputs = nuller_response(baseline,fov,sz,base_wavelength)
            pix2mas = fov*rad2mas/sz

            print("\nCalculating Exozodiacal")
            #exozodiacal flux (phot/s/m^2) per telescope
            exozodiacal = sf.calc_exozodiacal(star,outputs,local_exozodi,pix2mas,sz,spec)

            print("\nCalculating Leakage")
            #Calc stellar leakage flux (phot/s/m^2) per telescope
            leakage = stellar_leakage(star,nuller_response,baseline,base_wavelength)

            #pix2mas conversion, removing wavelength dependence
            #multiply by wavelength to get conversion factor for that wavelength
            wave_pix2mas = pix2mas/base_wavelength

            print("\nCalculating Signal")
            #signal flux (phot/s/m^2) per telescope
            signal = sf.calc_planet_signal(outputs,planet,wave_pix2mas,spec,mode)

            shot_noise = sf.calc_shot_noise(outputs,planet,wave_pix2mas,spec,mode)

            row_data = {"star_name":star.Name, "planet_name":planet.Name,
                        "universe_no":planet.UNumber,"star_no":star.SNumber,"planet_no":planet.PNumber,
                        "star_type":star.Stype,"star_distance (pc)":star.Dist,"baseline (m)":baseline,
                        "array_angle (mas)":star.PAngSep, "planet_angle (mas)":planet.PAngSep,
                        "star_flux (ph/s/m2)":star.flux,"planet_flux (ph/s/m2)":planet.flux,
                        "planet_temp (K)":planet.PTemp,"planet_radius (Earth_Rad)":planet.PRad,
                        "signal (ph/s/m2)":signal,
                        "shot (ph/s/m2)":shot_noise,
                        "leakage (ph/s/m2)":leakage,
                        "exozodiacal (ph/s/m2)":exozodiacal,
                        "zodiacal (ph/s)":zodiacal}

            ls_row_data.append(row_data)

    return ls_row_data
