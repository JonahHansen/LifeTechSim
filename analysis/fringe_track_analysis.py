#Quick analysis for calculating the fringe tracking rms needed to be photon limited

import numpy as np
from snr_calculator import load_results,grab_SNR_per_kernel,total_SNR
import json

#Limits for baselines
baseline_min = 5
baseline_max = 600

#Telescope parameters
D = 2
t = 3600
eta = 0.05

SNR_threshold = 7
n_universes = 10

#Output parameters
prefix = "/data/motley/jhansen/LifeSimData/avatar_run"
img_folder = "paper_plots/"
mode_names = ["Search", "Characterisation"]
arch_ls = [1,3,4,8,7,10]
n_scopes = [4,3,4,5,5,5]
arch_names = ["X-array","Kernel-3","Kernel-4","Kernel-5\n(0.66)","Kernel-5\n(1.03)","Kernel-5\n(1.68)"]


def calc_fringe_RMS(wave):

    overall_zod_ls = []
    fringe_rms_ls = []
    overall_wave_ls = []

    for a,n in zip(arch_ls,n_scopes):

        if a == 4:
            zod_fac = 0.5
        else:
            zod_fac = 1

        filename = prefix+"_"+str(a)+"_"+str(1)+"_"+str(wave)+".json"
        fin = open(filename,"r")
        data = json.load(fin)
        fin.close()

        waves = np.array(data["channel_centres (microns)"])
        print(waves)
        results = data["results"]

        zod_ls = []
        wave_ls = []
        temp_fringe_rms_ls = []

        for item in results:
            if (item['baseline'] <= baseline_max) &  (item['baseline'] >= baseline_min):
                snr_ls = grab_SNR_per_kernel(item,D,t,eta,zod_fac,True,n)
                tot_SNR = total_SNR(snr_ls)
                if tot_SNR > SNR_threshold:
                    idx = np.argmax(snr_ls)
                    A = np.pi*(D*np.sqrt(4/n))**2/4
                    zod = np.array(item["zodiacal"])[idx%50]/A/(np.array(item["leakage"]).flatten()[idx])
                    fringe_rms = np.sqrt(zod)*waves[idx%50]/2/np.pi
                    wave_ls.append(waves[idx%50])
                    zod_ls.append(zod)
                    temp_fringe_rms_ls.append(fringe_rms)

        overall_zod_ls.append(np.array(zod_ls))
        fringe_rms_ls.append(np.array(temp_fringe_rms_ls))
        overall_wave_ls.append(np.array(wave_ls))

    return overall_zod_ls,fringe_rms_ls,overall_wave_ls
