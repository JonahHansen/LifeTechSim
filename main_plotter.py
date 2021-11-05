import numpy as np
from snr_calculator import total_SNR_from_dict,load_results
import matplotlib.pyplot as plt
import json
import operator
from pokemon_matplotlib import pokemon_colours

baseline_min = 10
baseline_max = 600

D = 2
t = 3600

SNR_threshold = 7
n_universes = 10

prefix = "data/avatar_run"

arch = [1,3,4,7,8,9,10]
n_scopes = [4,3,4,5,5,5,5]
arch_names = ["Bracewell","Kernel 3","Kernel 4","Kernel 5 (1.03)","Kernel 5 (0.66)","Kernel 5 (2.67)","Kernel 5 (1.68)"]

#pokemon_colours("charmander")

def bar_plots(wave,eta):

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
    for a,n in zip(arch,n_scopes):

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

    plt.figure(2)
    plt.xticks(range(len(n_hab_ls)), arch_names)
    plt.xlabel('Architecture')
    plt.ylabel('Count (per universe)')
    plt.title('Habitable planets detected')
    plt.bar(range(len(n_hab_ls)), n_hab_ls)
    plt.legend()

    plt.figure(3)
    width = 0.3
    plt.xticks(range(len(n_rock_ls)), arch_names)
    plt.xlabel('Architecture')
    plt.ylabel('Count (per universe)')
    plt.title('Planets detected against radius/composition')
    plt.bar(np.arange(len(n_rock_ls))-width/2, n_rock_ls, width=width, label="Rocky (R<1.9)")
    plt.bar(np.arange(len(n_gas_ls))+ width/2, n_gas_ls, width=width, label="Gas (R>1.9)")
    plt.legend()

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

    plt.show()

    return

def n_max_arg(arr,n):
    ind = np.argpartition(arr,-n)[-n:]
    return ind[np.argsort(arr[ind])][::-1]

def sorted_indices(arr):
    a = [operator.itemgetter(0)(t) for t in sorted(enumerate(arr,1), key=operator.itemgetter(1))]
    a.reverse()
    return np.array(a)-1

def char_plots(wave,eta,n_planets):

    n_tot_ls = []

    bracewell_results = load_results(prefix,1,2,wave)
    bracewell_results.sort(key=lambda item: item.get("planet_name"))

    bracewell_results = [d for d in bracewell_results if d["habitable"] == "True"]

    bracewell_SNR_arr = []
    baseline_arr = []
    for item in bracewell_results:
        baseline_arr.append(item['baseline'])
        bracewell_SNR_arr.append(total_SNR_from_dict(item,D,t,eta,1,True,4))

    bracewell_SNR_arr = np.array(bracewell_SNR_arr)

    indices = sorted_indices(bracewell_SNR_arr)

    n_indices = []
    count = 0
    i = 0
    while count < n_planets:
        if (baseline_arr[indices[i]] <= baseline_max) and (baseline_arr[indices[i]] >= baseline_min):
            n_indices.append(indices[i])
            count += 1
        i += 1

    output_snr_ratio = []

    for a,n in zip(arch,n_scopes):

        if a == 4:
            zod_fac = 0.5
        else:
            zod_fac = 1

        results = load_results(prefix,a,2,wave)
        results.sort(key=lambda item: item.get("planet_name"))

        results = [d for d in results if d["habitable"] == "True"]

        ratio_SNR_arr = []

        for item,bracewell_SNR in zip(np.array(results)[n_indices],bracewell_SNR_arr[n_indices]):
            tot_SNR_arr = total_SNR_from_dict(item,D,t,eta,zod_fac,True,n)
            if (item['baseline'] > baseline_max) or (item['baseline'] < baseline_min):
                tot_SNR_arr = 0
            ratio = tot_SNR_arr/bracewell_SNR
            ratio_SNR_arr.append(ratio)

        output_snr_ratio.append(ratio_SNR_arr)

    plt.figure(1)
    plt.clf()
    #plt.xticks(range(len(output_snr_ratio[0])), arch_names)
    plt.xlabel('Planet no.')
    plt.ylabel('Relative SNR to Bracewell/Emma X-array')
    plt.title('Architecture relative SNR for characterisation of\n 30 highest SNR planets in Emma X-array configuration')
    for i,arch_ratio in enumerate(output_snr_ratio):
        plt.plot(np.array(range(n_planets))+1, arch_ratio, ls="",marker="_", mew=5,ms=20,label=arch_names[i])
    plt.legend()

    plt.figure(2)
    plt.clf()
    #plt.xticks(range(len(output_snr_ratio[0])), arch_names)
    plt.xlabel('Planet no.')
    plt.ylabel('Bracewell/Emma X-array SNR')
    plt.title('SNR for characterisation of\n 30 highest SNR planets in Emma X-array configuration')
    plt.plot(np.array(range(n_planets))+1, bracewell_SNR_arr[n_indices], ls="",marker="_", mew=5,ms=20)

    plt.show()
