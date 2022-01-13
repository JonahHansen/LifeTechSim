import numpy as np
from snr_calculator import total_SNR_from_dict,load_results,grab_SNR_per_kernel
import matplotlib.pyplot as plt
import matplotlib as mpl
import json
import operator
import cmasher as cmr
from collections import Counter

#List of colours
colours = cmr.take_cmap_colors('cmr.chroma', 6, cmap_range=(0.1,0.8), return_fmt='hex')

#Baseline range
baseline_min = 5
baseline_max = 600

#Telescope parameters
D = 2
t = 60*60*5
eta = 0.05

SNR_threshold = 7 #Threshold to determine if detectable
n_universes = 10

#Output parameters
prefix = "/data/motley/jhansen/LifeSimData/avatar_run"
img_folder = "paper_plots/"
mode_names = ["Search", "Characterisation"]
arch_ls = [1,3,4,8,7,10]
n_scopes = [4,3,4,5,5,5]
arch_names = ["X-array","Kernel-3","Kernel-4","Kernel-5\n(0.66)","Kernel-5\n(1.03)","Kernel-5\n(1.68)"]


"""
Plot the number of detected planets as a function of architecture (and other parameters) for a given reference wavelength
Inputs:
    wave = reference wavelength [10,15,18] (microns)
"""
def bar_plots(wave):

    n_tot_ls = [] #Total planets
    n_hab_ls = [] #Habitable zone planets
    n_hab_rock_ls = [] #Habitable, rocky and temperate planets

    planet_rad_divider = 1.5 #Separator between rocky and gaseous planets
    n_rock_ls = [] #Rocky planets
    n_gas_ls = [] #Gaseous planets

    temp_zone_min = 250 #Minimum temperate temperature
    temp_zone_max = 350 #Maximum temperate temperature
    n_cold_ls = [] #Cold planets
    n_temp_ls = [] #Temperate planets
    n_hot_ls = [] #Hot planets


    u_tot_ls = [] #Total planets
    u_hab_ls = [] #Habitable zone planets
    u_hab_rock_ls = [] #Habitable, rocky and temperate planets

    u_rock_ls = [] #Rocky planets
    u_gas_ls = [] #Gaseous planets

    u_cold_ls = [] #Cold planets
    u_temp_ls = [] #Temperate planets
    u_hot_ls = [] #Hot planets


    #For each architecture
    for a,n in zip(arch_ls,n_scopes):

        if a == 4:
            zod_fac = 0.5
        else:
            zod_fac = 1

        results = load_results(prefix,a,1,wave)

        tot_SNR_arr = []
        hab_SNR_arr = []
        hab_rock_SNR_arr = []

        rock_SNR_arr = []
        gas_SNR_arr = []

        cold_SNR_arr = []
        temp_SNR_arr = []
        hot_SNR_arr = []

        for item in results:
            if (item['baseline'] <= baseline_max) &  (item['baseline'] >= baseline_min):
                tot_SNR_arr.append((total_SNR_from_dict(item,D,t,eta,zod_fac,True,n),item["universe_no"]))
                if item["habitable"] == 'True':
                    hab_SNR_arr.append((total_SNR_from_dict(item,D,t,eta,zod_fac,True,n),item["universe_no"]))
                    if item["planet_radius"] <planet_rad_divider:
                        if item["planet_temp"] < temp_zone_max:
                            if item["planet_temp"] > temp_zone_min:
                                hab_rock_SNR_arr.append((total_SNR_from_dict(item,D,t,eta,zod_fac,True,n),item["universe_no"]))

                #Detections as a function of radii
                if item["planet_radius"] >planet_rad_divider:
                    gas_SNR_arr.append((total_SNR_from_dict(item,D,t,eta,zod_fac,True,n),item["universe_no"]))
                else:
                    rock_SNR_arr.append((total_SNR_from_dict(item,D,t,eta,zod_fac,True,n),item["universe_no"]))

                #Detections as a function of temperature
                if item["planet_temp"] > temp_zone_max:
                    hot_SNR_arr.append((total_SNR_from_dict(item,D,t,eta,zod_fac,True,n),item["universe_no"]))
                elif item["planet_temp"] < temp_zone_min:
                    cold_SNR_arr.append((total_SNR_from_dict(item,D,t,eta,zod_fac,True,n),item["universe_no"]))
                else:
                    temp_SNR_arr.append((total_SNR_from_dict(item,D,t,eta,zod_fac,True,n),item["universe_no"]))


        #Check if detected!
        temp_ls = [item[1] for item in tot_SNR_arr if item[0] > SNR_threshold]
        u_tot_ls.append(np.std(list(Counter(temp_ls).values())))
        n_tot_ls.append(np.sum(np.array([1 for item in tot_SNR_arr if item[0] > SNR_threshold]))/n_universes)

        temp_ls = [item[1] for item in hab_SNR_arr if item[0] > SNR_threshold]
        u_hab_ls.append(np.std(list(Counter(temp_ls).values())))
        n_hab_ls.append(np.sum(np.array([1 for item in hab_SNR_arr if item[0] > SNR_threshold]))/n_universes)

        temp_ls = [item[1] for item in hab_rock_SNR_arr if item[0] > SNR_threshold]
        u_hab_rock_ls.append(np.std(list(Counter(temp_ls).values())))
        n_hab_rock_ls.append(np.sum(np.array([1 for item in hab_rock_SNR_arr if item[0] > SNR_threshold]))/n_universes)

        temp_ls = [item[1] for item in rock_SNR_arr if item[0] > SNR_threshold]
        u_rock_ls.append(np.std(list(Counter(temp_ls).values())))
        n_rock_ls.append(np.sum(np.array([1 for item in rock_SNR_arr if item[0] > SNR_threshold]))/n_universes)

        temp_ls = [item[1] for item in gas_SNR_arr if item[0] > SNR_threshold]
        u_gas_ls.append(np.std(list(Counter(temp_ls).values())))
        n_gas_ls.append(np.sum(np.array([1 for item in gas_SNR_arr if item[0] > SNR_threshold]))/n_universes)

        temp_ls = [item[1] for item in cold_SNR_arr if item[0] > SNR_threshold]
        u_cold_ls.append(np.std(list(Counter(temp_ls).values())))
        n_cold_ls.append(np.sum(np.array([1 for item in cold_SNR_arr if item[0] > SNR_threshold]))/n_universes)

        temp_ls = [item[1] for item in temp_SNR_arr if item[0] > SNR_threshold]
        u_temp_ls.append(np.std(list(Counter(temp_ls).values())))
        n_temp_ls.append(np.sum(np.array([1 for item in temp_SNR_arr if item[0] > SNR_threshold]))/n_universes)

        temp_ls = [item[1] for item in hot_SNR_arr if item[0] > SNR_threshold]
        u_hot_ls.append(np.std(list(Counter(temp_ls).values())))
        n_hot_ls.append(np.sum(np.array([1 for item in hot_SNR_arr if item[0] > SNR_threshold]))/n_universes)

        """
        n_tot_ls.append(np.sum(np.array(tot_SNR_arr)>SNR_threshold)/n_universes)
        n_hab_ls.append(np.sum(np.array(hab_SNR_arr)>SNR_threshold)/n_universes)
        n_hab_rock_ls.append(np.sum(np.array(hab_rock_SNR_arr)>SNR_threshold)/n_universes)

        n_rock_ls.append(np.sum(np.array(rock_SNR_arr)>SNR_threshold)/n_universes)
        n_gas_ls.append(np.sum(np.array(gas_SNR_arr)>SNR_threshold)/n_universes)

        n_cold_ls.append(np.sum(np.array(cold_SNR_arr)>SNR_threshold)/n_universes)
        n_temp_ls.append(np.sum(np.array(temp_SNR_arr)>SNR_threshold)/n_universes)
        n_hot_ls.append(np.sum(np.array(hot_SNR_arr)>SNR_threshold)/n_universes)
        """

    print(u_tot_ls)

    plt.figure(1)
    plt.clf()
    plt.xticks(range(len(n_tot_ls)), arch_names)
    plt.xlabel('Architecture')
    plt.ylabel('Average Count')
    plt.title('Total planets detected')
    plt.bar(range(len(n_tot_ls)), n_tot_ls, color=colours[0], yerr=u_tot_ls, ecolor="#cccccc")
    plt.legend()
    plt.savefig(img_folder+"Total_planets_bar_%s_micron.pdf"%wave,bbox_inches='tight',dpi=100)

    plt.figure(2)
    plt.clf()
    plt.xticks(range(len(n_hab_ls)), arch_names)
    plt.xlabel('Architecture')
    plt.ylabel('Average Count')
    plt.title('Planets detected in the habitable zone')
    plt.bar(range(len(n_hab_ls)), n_hab_ls, color=colours[0], yerr=u_hab_ls, ecolor="#cccccc")
    plt.legend()
    plt.savefig(img_folder+"Habitable_planets_bar_%s_micron.pdf"%wave,bbox_inches='tight',dpi=100)

    plt.figure(3)
    plt.clf()
    width = 0.3
    plt.xticks(range(len(n_rock_ls)), arch_names)
    plt.xlabel('Architecture')
    plt.ylabel('Average Count')
    plt.title('Planets detected against radius/composition')
    plt.bar(np.arange(len(n_rock_ls))-width/2, n_rock_ls, width=width, label="Rocky (R_E<1.9)",color=colours[0], yerr=u_rock_ls, ecolor="#cccccc")
    plt.bar(np.arange(len(n_gas_ls))+ width/2, n_gas_ls, width=width, label="Gas (R_E>1.9)",color=colours[2], yerr=u_gas_ls, ecolor="#cccccc")
    plt.legend()
    plt.savefig(img_folder+"Radius_planets_bar_%s_micron.pdf"%wave,bbox_inches='tight',dpi=100)

    plt.figure(4)
    plt.clf()
    width = 0.25
    plt.xticks(range(len(n_cold_ls)), arch_names)
    plt.xlabel('Architecture')
    plt.ylabel('Average Count')
    plt.title('Planets detected against temperature')
    plt.bar(np.arange(len(n_cold_ls)) - width, n_cold_ls, width=width, label="Cold (T<250K)",color=colours[0], yerr=u_cold_ls, ecolor="#cccccc")
    plt.bar(np.arange(len(n_temp_ls)), n_temp_ls, width=width, label="Temperate (250<T<350K)",color=colours[2], yerr=u_temp_ls, ecolor="#cccccc")
    plt.bar(np.arange(len(n_hot_ls))+ width, n_hot_ls, width=width, label="Hot (T>350K)",color=colours[4], yerr=u_hot_ls, ecolor="#cccccc")
    plt.legend()
    plt.savefig(img_folder+"Temperature_planets_bar_%s_micron.pdf"%wave,bbox_inches='tight',dpi=100)

    plt.figure(5)
    plt.clf()
    plt.xticks(range(len(n_hab_rock_ls)), arch_names)
    plt.xlabel('Architecture')
    plt.ylabel('Average Count')
    plt.title('Temperate rocky planets detected in the habitable zone')
    plt.bar(range(len(n_hab_rock_ls)), n_hab_rock_ls, color=colours[0], yerr=u_hab_rock_ls, ecolor="#cccccc")
    plt.legend()
    plt.savefig(img_folder+"Habitable_temperate_rocky_planets_bar_%s_micron.pdf"%wave,bbox_inches='tight',dpi=100)

    return


# Get a list of sorted indices based on an input array
def sorted_indices(arr):
    a = [operator.itemgetter(0)(t) for t in sorted(enumerate(arr,1), key=operator.itemgetter(1))]
    a.reverse()
    return np.array(a)-1


"""
Characterisation plot - plot of the top N habitable planets at a reference wavelength
and their relative SNRs against a given architecture for all other architectures

Inputs:
    wave = reference wavelength [10,15,18] (microns)
    base_arch = index of the architecture to calculate reference SNR
    n_planets = how many planets in plot?

"""
def char_plots(wave,base_arch,n_planets=25):

    n_tot_ls = []

    if base_arch == 4:
        zod_fac = 0.5
    else:
        zod_fac = 1

    #Number of telescopes
    if base_arch >= 7:
        n = 5
    elif base_arch == 3:
        n = 3
    else:
        n = 4

    base_results = load_results(prefix,base_arch,2,wave)
    base_results.sort(key=lambda item: item.get("planet_name"))

    base_results = [d for d in base_results if (d["habitable"] == "True")]

    base_SNR_arr = []
    baseline_arr = []
    for item in base_results:
        baseline_arr.append(item['baseline'])
        base_SNR_arr.append(total_SNR_from_dict(item,D,t,eta,zod_fac,True,n))

    base_SNR_arr = np.array(base_SNR_arr)

    #indices of the planets sorted by SNR using the base_architecture
    indices = sorted_indices(base_SNR_arr)

    #Take top n planets, and ensure that all are within the baseline limits
    n_indices = []
    count = 0
    i = 0
    while count < n_planets:
        if (baseline_arr[indices[i]] <= baseline_max) and (baseline_arr[indices[i]] >= baseline_min):
            n_indices.append(indices[i])
            count += 1
        i += 1

    output_snr_ratio = []

    #Calculate relative SNR for each architecture
    for a,n in zip(arch_ls,n_scopes):

        if a == 4:
            zod_fac = 0.5
        else:
            zod_fac = 1

        results = load_results(prefix,a,2,wave)
        results.sort(key=lambda item: item.get("planet_name"))

        results = [d for d in results if d["habitable"] == "True"]

        ratio_SNR_arr = []

        for item,base_SNR in zip(np.array(results)[n_indices],base_SNR_arr[n_indices]):
            tot_SNR_arr = total_SNR_from_dict(item,D,t,eta,zod_fac,True,n)
            if (item['baseline'] > baseline_max) or (item['baseline'] < baseline_min):
                tot_SNR_arr = 0
            ratio = tot_SNR_arr/base_SNR
            ratio_SNR_arr.append(ratio)

        output_snr_ratio.append(ratio_SNR_arr)

    plt.figure(1)
    plt.clf()
    #plt.xticks(range(len(output_snr_ratio[0])), arch_names)
    plt.xlabel('Planet no.')
    plt.ylabel('Relative SNR to X-array')
    plt.title('Architecture relative SNR for characterisation of \nthe 25 highest SNR planets in X-array configuration')
    for i,arch_ratio in enumerate(output_snr_ratio):
        plt.plot(np.array(range(n_planets))+1, arch_ratio, ls="",marker="_", mew=3,ms=10,label=arch_names[i],color=colours[i])
    plt.legend()
    plt.savefig(img_folder+"Char_plot_%s_arch_%s_micron.pdf"%(base_arch,wave),bbox_inches='tight',dpi=100)

    return


"""
Plot all architectures SNR as a function of wavelength for a given planet

Inputs:
    mode = search (0) or characterisation (1)
    wave = reference wavelength [10,15,18]
    planet_index = index of (habitable) planet to plot for

"""
def snr_wave_plot(mode,wave,planet_index):

    plt.figure(1)
    plt.clf()
    for i,a,n,name in zip(range(6),arch_ls,n_scopes,arch_names):

        if a == 4:
            zod_fac = 0.5
        else:
            zod_fac = 1

        results = load_results(prefix,a,mode,wave)
        results.sort(key=lambda item: item.get("planet_name"))

        results = [d for d in results if d["habitable"] == "True"]

        item = np.array(results)[planet_index]

        #print planet data
        print("###################\nPlanet data\n####################")
        print(item["planet_name"])
        print("Star distance = %s"%item["star_distance"])
        print("Star type = %s"%item["star_type"])
        print("Planet angle = %s"%item["planet_angle"])
        print("Array angle = %s"%item["array_angle"])
        print("Array baseline = %s"%item["baseline"])
        print("Planet temp = %s"%item["planet_temp"])
        print("Planet radius = %s"%item["planet_radius"])
        print("Habitable? = %s"%item["habitable"])

        snr = grab_SNR_per_kernel(item,D,t,eta,zod_fac,True,n)

        combined_snr = np.sqrt(np.sum(snr**2,axis=0))

        waves = np.linspace(4,19,50)

        plt.plot(waves,combined_snr,label=name,color=colours[i])

    plt.xlabel(r"Wavelength ($\mu$m)")
    plt.ylabel("SNR")
    plt.legend()
    plt.savefig(img_folder+"SNR_wavelength_1_18_%s.pdf"%planet_index,bbox_inches='tight',dpi=100)

    return

########################### TEST PLOTS ###########################################

"""
Plot characterisation SNR as a function of wavelength for a given planet with certain components removed

Inputs:
    arch = architecure index
    wave = wavelength [10,15,18] (microns)
    n_telescopes = number of telescopes used in given architecture
    planet_index = index of planet to plot for

"""
def snr_component_plot(arch,wave,n_telescopes,planet_index):

    output_snr_ratio = []

    if arch == 4:
        zod_fac = 0.5
    else:
        zod_fac = 1

    results = load_results(prefix,arch,2,wave)
    results.sort(key=lambda item: item.get("planet_name"))

    results = [d for d in results if d["habitable"] == "True"]

    item = np.array(results)[planet_index]

    print("###################\nPlanet data\n####################")
    print("Star distance = %s"%item["star_distance"])
    print("Star type = %s"%item["star_type"])
    print("Planet angle = %s"%item["planet_angle"])
    print("Array angle = %s"%item["array_angle"])
    print("Array baseline = %s"%item["baseline"])
    print("Planet temp = %s"%item["planet_temp"])
    print("Planet radius = %s"%item["planet_radius"])
    print("Habitable? = %s"%item["habitable"])

    #Default SNR
    snr_1 = grab_SNR_per_kernel(item,D,t,eta,zod_fac,True,n_telescopes)
    #No zodiacal
    snr_2 = grab_SNR_per_kernel(item,D,t,eta,0,True,n_telescopes)
    #No exozodiacal
    snr_3 = grab_SNR_per_kernel(item,D,t,eta,zod_fac,True,n_telescopes,exozodfac=0)
    #No stellar leakage
    snr_4 = grab_SNR_per_kernel(item,D,t,eta,zod_fac,True,n_telescopes,stellarfac=0)

    linestyles = ["-","--"]

    waves = np.linspace(4,19,50)
    plt.figure(1)
    plt.clf()
    for j in range(len(snr_1)):
        plt.plot(waves,snr_2[j],c="g",ls=linestyles[j],label="K%s, No zodi noise"%(j+1))
        plt.plot(waves,snr_3[j],c="r",ls=linestyles[j],label="K%s, No exozodi noise"%(j+1))
        plt.plot(waves,snr_4[j],c="c",ls=linestyles[j],label="K%s, No leakage noise"%(j+1))
        plt.plot(waves,snr_1[j],c="b",ls=linestyles[j],label="K%s, All noise"%(j+1))

    plt.xlabel("Wavelength (microns)")
    plt.ylabel("SNR")
    plt.legend()
    plt.show()
    return


"""
Helper function to round number to a given number of sig figs

Inputs:
    x - number
    p - precision (number of sig figs)
"""
def round_sig_figs(x, p):
    x_positive = np.where(np.isfinite(x) & (x != 0), np.abs(x), 10**(p-1))
    mags = 10 ** (p - 1 - np.floor(np.log10(x_positive)))
    return np.round(x * mags) / mags


"""
Plot noise sources for a given planet as a function of wavelength

Inputs:
    arch = architecure index
    mode = search (0) or characterisation (1)
    wave = wavelength [10,15,18] (microns)
    D = telescope diameter (m)
    n_telescopes = number of telescopes used in given architecture
    planet_index = index of planet to plot for

"""
def noise_contributions_plot(arch,mode,wave,D,n_telescopes,planet_index):
    results = load_results(prefix,arch,mode,wave)
    item = results[planet_index]

    print("###################\nPlanet data %s\n####################"%arch)
    print("Star distance = %s"%item["star_distance"])
    print("Star type = %s"%item["star_type"])
    print("Planet angle = %s"%item["planet_angle"])
    print("Array angle = %s"%item["array_angle"])
    print("Array baseline = %s"%item["baseline"])
    print("Planet temp = %s"%item["planet_temp"])
    print("Planet radius = %s"%item["planet_radius"])
    print("Habitable? = %s"%item["habitable"])

    D *= np.sqrt(4/n_telescopes)

    A = np.pi*D**2/4

    if arch == 4:
        zod_fac = 0.5
    else:
        zod_fac = 1

    linestyles = ["-","--",".."]

    signal = 2*np.array(item["signal"])*A
    shot = 2*np.array(item["shot"])*A
    leakage = 2*np.array(item["leakage"])*A
    exozodiacal = 2*np.array(item["exozodiacal"])*A
    zodiacal = 2*np.array(item["zodiacal"])*zod_fac

    waves = np.linspace(4,19,50)

    plt.figure(1)
    plt.clf()

    for j in range(len(signal)):
        plt.plot(waves,np.log10(signal[j]),c="k",ls=linestyles[j],label="K%s, Signal"%(j+1))
        plt.plot(waves,np.log10(shot[j]),c="r",ls=linestyles[j],label="K%s, Shot noise"%(j+1))
        plt.plot(waves,np.log10(leakage[j]),c="c",ls=linestyles[j],label="K%s, Stellar leakage noise"%(j+1))
        plt.plot(waves,np.log10(exozodiacal[j]),c="b",ls=linestyles[j],label="K%s, Exozodiacal noise"%(j+1))
        if j ==0:
            plt.plot(waves,np.log10(zodiacal),c="g",ls=linestyles[j],label="K%s, Zodiacal noise"%(j+1))

    plt.xlabel("Wavelength (um)")
    plt.ylabel("Photon rate (ph/s)")
    plt.title("Photon rates for the various signal and noise components\n for the %s architecture, %s mode\n reference wavelength of %s and D=%s"%(arch,mode,wave,round_sig_figs(D,3)))
    plt.legend()
    plt.show()
    return
