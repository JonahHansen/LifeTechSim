import numpy as np
from snr_calculator import total_SNR_from_dict,load_results,grab_SNR_per_kernel,get_data_one_planet
import matplotlib.pyplot as plt
import matplotlib as mpl
import json
import operator
import cycler
#from pokemon_matplotlib import pokemon_colours

baseline_min = 5 #FINE FOR ALL (As baseline refers to smallest baseline)
baseline_max = 600 #Possibly problematic... especially for bracewell which has 6x smaller limit

D = 2
t = 3600

SNR_threshold = 7
n_universes = 10

eta = 0.05

prefix = "/data/motley/jhansen/LifeSimData/avatar_run"
img_folder = "plots/"
mode_names = ["Search", "Characterisation"]
arch_ls = [1,2,3,4,7,8,9,10]
n_scopes = [4,4,3,4,5,5,5,5]
arch_names = ["Bracewell","Linear","Kernel 3","Kernel 4","Kernel 5 (1.03)","Kernel 5 (0.66)","Kernel 5 (2.67)","Kernel 5 (1.68)"]

#pokemon_colours("charmander")

def bar_plots(wave):

    n_tot_ls = []
    n_hab_ls = []

    planet_rad_divider = 1.9
    n_rock_ls = []
    n_gas_ls = []

    temp_zone_min = 250
    temp_zone_max = 350
    n_cold_ls = []
    n_temp_ls = []
    n_hot_ls = []
    for a,n in zip(arch_ls,n_scopes):

        if a == 4:
            zod_fac = 0.5
        else:
            zod_fac = 1

        results = load_results(prefix,a,1,wave)

        tot_SNR_arr = []
        hab_SNR_arr = []

        rock_SNR_arr = []
        gas_SNR_arr = []

        cold_SNR_arr = []
        temp_SNR_arr = []
        hot_SNR_arr = []

        for item in results:
            if (item['baseline'] <= baseline_max) &  (item['baseline'] >= baseline_min):
                tot_SNR_arr.append(total_SNR_from_dict(item,D,t,eta,zod_fac,True,n))
                if item["habitable"] == 'True':
                    hab_SNR_arr.append(total_SNR_from_dict(item,D,t,eta,zod_fac,True,n))

                if item["planet_radius"] >planet_rad_divider:
                    gas_SNR_arr.append(total_SNR_from_dict(item,D,t,eta,zod_fac,True,n))
                else:
                    rock_SNR_arr.append(total_SNR_from_dict(item,D,t,eta,zod_fac,True,n))

                if item["planet_temp"] > temp_zone_max:
                    hot_SNR_arr.append(total_SNR_from_dict(item,D,t,eta,zod_fac,True,n))
                elif item["planet_temp"] < temp_zone_min:
                    cold_SNR_arr.append(total_SNR_from_dict(item,D,t,eta,zod_fac,True,n))
                else:
                    temp_SNR_arr.append(total_SNR_from_dict(item,D,t,eta,zod_fac,True,n))


        n_tot_ls.append(np.sum(np.array(tot_SNR_arr)>SNR_threshold)/n_universes)
        n_hab_ls.append(np.sum(np.array(hab_SNR_arr)>SNR_threshold)/n_universes)

        n_rock_ls.append(np.sum(np.array(rock_SNR_arr)>SNR_threshold)/n_universes)
        n_gas_ls.append(np.sum(np.array(gas_SNR_arr)>SNR_threshold)/n_universes)

        n_cold_ls.append(np.sum(np.array(cold_SNR_arr)>SNR_threshold)/n_universes)
        n_temp_ls.append(np.sum(np.array(temp_SNR_arr)>SNR_threshold)/n_universes)
        n_hot_ls.append(np.sum(np.array(hot_SNR_arr)>SNR_threshold)/n_universes)

    plt.figure(1)
    plt.xticks(range(len(n_tot_ls)), arch_names)
    plt.xlabel('Architecture')
    plt.ylabel('Count (per universe)')
    plt.title('Total planets detected')
    plt.bar(range(len(n_tot_ls)), n_tot_ls)
    plt.legend()
    plt.savefig(img_folder+"Total_planets_bar_%s_micron.png"%wave)

    plt.figure(2)
    plt.xticks(range(len(n_hab_ls)), arch_names)
    plt.xlabel('Architecture')
    plt.ylabel('Count (per universe)')
    plt.title('Habitable planets detected')
    plt.bar(range(len(n_hab_ls)), n_hab_ls)
    plt.legend()
    plt.savefig(img_folder+"Habitable_planets_bar_%s_micron.png"%wave)

    plt.figure(3)
    width = 0.3
    plt.xticks(range(len(n_rock_ls)), arch_names)
    plt.xlabel('Architecture')
    plt.ylabel('Count (per universe)')
    plt.title('Planets detected against radius/composition')
    plt.bar(np.arange(len(n_rock_ls))-width/2, n_rock_ls, width=width, label="Rocky (R<1.9)")
    plt.bar(np.arange(len(n_gas_ls))+ width/2, n_gas_ls, width=width, label="Gas (R>1.9)")
    plt.legend()
    plt.savefig(img_folder+"Radius_planets_bar_%s_micron.png"%wave)

    plt.figure(4)
    width = 0.25
    plt.xticks(range(len(n_cold_ls)), arch_names)
    plt.xlabel('Architecture')
    plt.ylabel('Count (per universe)')
    plt.title('Planets detected against temperature')
    plt.bar(np.arange(len(n_cold_ls)) - width, n_cold_ls, width=width, label="Cold (T<250K)")
    plt.bar(np.arange(len(n_temp_ls)), n_temp_ls, width=width, label="Temperate (250<T<350K)")
    plt.bar(np.arange(len(n_hot_ls))+ width, n_hot_ls, width=width, label="Hot (T>350K)")
    plt.legend()
    plt.savefig(img_folder+"Temperature_planets_bar_%s_micron.png"%wave)

    #plt.show()

    return

def sorted_indices(arr):
    a = [operator.itemgetter(0)(t) for t in sorted(enumerate(arr,1), key=operator.itemgetter(1))]
    a.reverse()
    return np.array(a)-1

def char_plots(wave,base_arch,n_planets):

    n_tot_ls = []

    base_results = load_results(prefix,base_arch,2,wave)
    base_results.sort(key=lambda item: item.get("planet_name"))

    #bracewell_results = [d for d in bracewell_results if (d["habitable"] == "True") & (d["planet_radius"] < 1.9)]
    base_results = [d for d in base_results if (d["habitable"] == "True")]

    base_SNR_arr = []
    baseline_arr = []
    for item in base_results:
        baseline_arr.append(item['baseline'])
        base_SNR_arr.append(total_SNR_from_dict(item,D,t,eta,1,True,4))

    base_SNR_arr = np.array(base_SNR_arr)

    indices = sorted_indices(base_SNR_arr)

    n_indices = []
    count = 0
    i = 0
    while count < n_planets:
        if (baseline_arr[indices[i]] <= baseline_max) and (baseline_arr[indices[i]] >= baseline_min):
            n_indices.append(indices[i])
            count += 1
        i += 1

    output_snr_ratio = []

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
    plt.ylabel('Relative SNR to Bracewell/Emma X-array')
    plt.title('Architecture relative SNR for characterisation at %s um of\n %s highest SNR planets in Emma X-array configuration'%(wave,n_planets))
    for i,arch_ratio in enumerate(output_snr_ratio):
        plt.plot(np.array(range(n_planets))+1, arch_ratio, ls="",marker="_", mew=5,ms=20,label=arch_names[i])
    plt.legend()
    plt.savefig("Char_plot_%s_arch_%s_micron.png"%(base_arch,wave))

    """
    plt.figure(2)
    plt.clf()
    #plt.xticks(range(len(output_snr_ratio[0])), arch_names)
    plt.xlabel('Planet no.')
    plt.ylabel('Bracewell/Emma X-array SNR')
    plt.title('SNR for characterisation at %s um of\n %s highest SNR planets in Emma X-array configuration'%(wave,n_planets))
    plt.plot(np.array(range(n_planets))+1, base_SNR_arr[n_indices], ls="",marker="_", mew=5,ms=20)
    """
    #plt.show()
    return

#Plot SNR as a function of wavelength for a given planet with certain components removed
def snr_component_plot(arch,n_telescopes,wave,planet_index):

    output_snr_ratio = []

    if arch == 4:
        zod_fac = 0.5
    else:
        zod_fac = 1

    results = load_results(prefix,arch,2,wave)
    results.sort(key=lambda item: item.get("planet_name"))

    results = [d for d in results if d["habitable"] == "True"]

    #plt.ioff()

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

    snr_1 = grab_SNR_per_kernel(item,D,t,eta,zod_fac,True,n_telescopes)
    snr_2 = grab_SNR_per_kernel(item,D,t,eta,0,True,n_telescopes)
    snr_3 = grab_SNR_per_kernel(item,D,t,eta,zod_fac,True,n_telescopes,exozodfac=0)
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

    plt.legend()
    plt.show()
    return

#Plot all architectures SNR as a function of wavelength
def snr_wave_plot(mode,wave,planet_index):

    color = plt.cm.tab10(np.linspace(0, 1,10))
    mpl.rcParams['axes.prop_cycle'] = cycler.cycler('color', color)

    plt.figure(1)
    plt.clf()
    for a,n,name in zip(arch_ls,n_scopes,arch_names):

        if a == 4:
            zod_fac = 0.5
        else:
            zod_fac = 1

        results = load_results(prefix,a,mode,wave)
        results.sort(key=lambda item: item.get("planet_name"))

        results = [d for d in results if d["habitable"] == "True"]

        #plt.ioff()

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


        snr = grab_SNR_per_kernel(item,D,t,eta,zod_fac,True,n)

        combined_snr = np.sqrt(np.sum(snr**2,axis=0))

        waves = np.linspace(4,19,50)

        plt.plot(waves,combined_snr,label=name)

    plt.xlabel("Wavelength (um)")
    plt.ylabel("SNR")
    plt.title("SNR against wavelength for %s\n with a reference wavelength of %s um"%(mode_names[mode-1],wave))
    plt.legend()
    plt.show()
    return


def round_sig_figs(x, p):
    x_positive = np.where(np.isfinite(x) & (x != 0), np.abs(x), 10**(p-1))
    mags = 10 ** (p - 1 - np.floor(np.log10(x_positive)))
    return np.round(x * mags) / mags


#Plot noise sources for a given planet as a function of wavelength
def noise_contributions_plot(planet_index,arch,mode,wave,D,num_telescopes):
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

    D *= np.sqrt(4/num_telescopes)

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

def make_plots(wave):
    bar_plots(wave)
    char_plots(wave,1,25)
    char_plots(wave,2,25)
    char_plots(wave,3,25)
    char_plots(wave,4,25)
    char_plots(wave,7,25)
    char_plots(wave,8,25)
    char_plots(wave,9,25)
    char_plots(wave,10,25)
