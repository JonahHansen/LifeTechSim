#File to plot the outputs from the dist_dependence_sim simulation.py

import numpy as np
from snr_calculator import total_SNR_from_dict
import matplotlib.pyplot as plt
import json
from bisect import bisect_left, bisect_right
import cmasher as cmr

star_types = ["A6V","F7V","G2V","K2V","M5V"] #List of star typesa
planets = ["Inner HZ", "Mid HZ", "Outer HZ"] #List of planet locations

colours = cmr.take_cmap_colors('cmr.chroma', 6, cmap_range=(0.1,0.8), return_fmt='hex')

#Min and max baselines allowed
baseline_min = 5
baseline_max = 600

#Telescope parameters
D = 2
t = 3600
eta = 0.05

SNR_threshold = 7

data_prefix = "/data/motley/jhansen/LifeSimData/avatar_dist_folder/avatar_dist_run"

img_folder = "paper_plots/"

arch_ls = [1,3,4,8,7,10] #indices of architecture run
n_scopes = [4,3,4,5,5,5] #number of telescopes for each architecture
arch_names = ["X-array","Kernel-3","Kernel-4","Kernel-5 (0.66)","Kernel-5 (1.03)","Kernel-5 (1.68)"]

"""
Load list of dictionaries for a given star and planet from a JSON file.
Indices are based on how the output files are named.


Inputs:
    prefix = prefix of JSON files
    arch = architecure index
    wave = wavelength [10,15,18] (microns)
    star_index = index of "A6V","F7V","G2V","K2V","M5V" [1,2,3,4,5]
    planets = index of "InnerHZ", "MidHZ", "OutHZ" [1,2,3]

Output:
    list of dictionaries
"""
def load_results(arch,wave,star_index,planet_index):

    filename = data_prefix+"_"+str(arch)+"_"+str(wave)+".json"
    fin = open(filename,"r")
    data = json.load(fin)
    fin.close()

    results = data["results"]

    def condition(dict):
        return (dict["star_type"] == star_types[star_index]) & (dict["planet_no"] == planet_index+1)

    filtered_list = [d for d in results if condition(d)]

    sorted_list = sorted(filtered_list, key=lambda d: d['star_distance'])

    return sorted_list

"""
Get a list of SNRs, along with their baselines, for each distance of a given star and planet

Inputs:
    arch = architecure index
    wave = wavelength [10,15,18] (microns)
    star_index = index of "A6V","F7V","G2V","K2V","M5V" [1,2,3,4,5]
    planets = index of "InnerHZ", "MidHZ", "OutHZ" [1,2,3]
    n_telescopes = number of telescopes used in given architecture

Outputs:
    list of distances in pc
    list of baselines in m
    list of SNRs

"""
def get_snr_by_dist(arch,wave,star_index,planet_index,n_telescopes):
    #load results
    res = load_results(arch,wave,star_index,planet_index)

    dist = []
    baseline = []
    snr = []
    for d in res:
        baseline.append(d['baseline'])
        dist.append(d['star_distance'])
        if arch == 4:
            snr.append(total_SNR_from_dict(d,D,t,eta,0.5,const_total_area=True,num_telescopes=n_telescopes))
        else:
            snr.append(total_SNR_from_dict(d,D,t,eta,const_total_area=True,num_telescopes=n_telescopes))

    return np.array(dist), np.array(baseline), np.array(snr)


"""
Main function to plot the SNR data as a function of architecture, both in absolute and relative SNRs against the X-array configuration.

Inputs:
    wave = wavelength [10,15,18] (microns)
    star_index = index of "A6V","F7V","G2V","K2V","M5V" [1,2,3,4,5]
    planets = index of "InnerHZ", "MidHZ", "OutHZ" [1,2,3]

Outputs:
    Plot of SNR against distance for each architecture
    Plot of relative SNR (to X-array) against distance for each architecture

"""
def make_plot(wave,star_index,planet_index):

    #Calculate x-array (bracewell) SNR
    dist_a,baseline,snr_b = get_snr_by_dist(1,wave,star_index,planet_index,4)

    plt.figure(1)
    plt.clf()
    plt.figure(2)
    plt.clf()

    #Plot each architecture's data, making baselines outside of set range as a dashed line
    for a,n,name,c in zip(arch,n_scopes,arch_names,colours):
        dist, baseline, snr = get_snr_by_dist(a,wave,star_index,planet_index,n)
        bad_index_max = bisect_right(baseline-baseline_max,0)
        bad_index_min = bisect_left(baseline-baseline_min,0)
        plt.figure(1)
        plt.plot(dist[:bad_index_min],np.log10(snr[:bad_index_min]),label=name,linestyle="--",color=c)
        plt.plot(dist[bad_index_max:],np.log10(snr[bad_index_max:]),label="",linestyle="--",color=c)
        plt.plot(dist[bad_index_min:bad_index_max],np.log10(snr[bad_index_min:bad_index_max]),label="",linestyle="-",color=c)

        plt.figure(2)
        plt.plot(dist[:bad_index_min],snr[:bad_index_min]/snr_b[:bad_index_min],label=name,linestyle="--",color=c)
        plt.plot(dist[bad_index_max:],snr[bad_index_max:]/snr_b[bad_index_max:],label="",linestyle="--",color=c)
        plt.plot(dist[bad_index_min:bad_index_max],snr[bad_index_min:bad_index_max]/snr_b[bad_index_min:bad_index_max],label="",linestyle="-",color=c)


    plt.figure(1)
    plt.axhline(y = np.log10(SNR_threshold), color = 'k', linestyle = '--')
    plt.xlabel("Distance (pc)")
    plt.ylabel("[SNR]")
    plt.title("SNR against distance for planets in the %s\n for a %s type star at %s microns" %(planets[planet_index],star_types[star_index],wave))
    plt.legend()

    plt.figure(2)
    plt.xlabel("Distance (pc)")
    plt.ylabel("Relative SNR to X-array")
    plt.title("Relative SNR against distance for planets in the %s\n for a %s type star" %(planets[planet_index],star_types[star_index]))
    plt.legend()
    plt.savefig(img_folder+"Dist_plot_%s_planet_%s_star.pdf"%(planets[planet_index],star_types[star_index]),bbox_inches='tight',dpi=100)

    return
