import json
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
plt.ion()


"""
Calculate the signal to noise from the various photon amounts

Inputs:
    signal, shot, leakage, zodiacal, exozodiacal = number of photons of each contribution

Outputs:
    SNR
"""
def SNR(signal,shot,leakage,zodiacal,exozodiacal):
    sig = signal
    noise = 2*shot + 2*leakage + 2*zodiacal + 2*exozodiacal #Factors of two => difference of two responses
    return sig/np.sqrt(noise)


"""
Calculate the signal to noise of one planet, for each kernel and wavelength channel

Inputs:
    dict = input result dictionary of a planet (from the main simulation)
    D = telescope diameter (m)
    t = exposure time (s)
    eta = throughput

Outputs:
    List of SNRs for each kernel and wavelength
"""
def grab_SNR_per_kernel(dict,D,t,eta):

    A = np.pi*D**2/4

    signal = np.array(dict["signal (ph/s/m2)"])*A*t*eta
    shot = np.array(dict["shot (ph/s/m2)"])*A*t*eta
    leakage = np.array(dict["leakage (ph/s/m2)"])*A*t*eta
    exozodiacal = np.array(dict["exozodiacal (ph/s/m2)"])*A*t*eta
    zodiacal = np.array(dict["zodiacal (ph/s)"])*t*eta

    return SNR(signal, shot, leakage, zodiacal, exozodiacal)


"""
Take a list of all the SNRs for each kernel and wavelength and combine them

Input:
    SNR_array = list of SNRs

Output:
    Total SNR
"""
def total_SNR(SNR_array):
    return np.sqrt(np.sum(SNR_array**2))


"""
Rough cost of the telescope array from Stahl 2020

Inputs:
    N = number of telescopes
    D = diameter of telescopes (m)
    lam = diffraction limited wavelength (um)
    T = temperature of the array (K)

Output:
    Cost in USD millions
"""
def cost(N,D,lam,T=20):
    return N*300*D**1.7*lam**(-0.5)*T**(-0.25)


"""
Load list of dictionaries from a JSON file.
Indices are based on how the output files are named

Inputs:
    prefix = prefix of JSON files
    arch = architecure index [1,2,3,4,5,6,7,8]
    mode = mode index [1,2]
    wave = wavelength [10,15,18]

Output:
    list of dictionaries
"""
def load_results(prefix,arch,mode,wave):
    filename = prefix+"_"+str(arch)+"_"+str(mode)+"_"+str(wave)+".json"
    fin = open(filename,"r")
    data = json.load(fin)
    fin.close()
    return data["results"]


"""
Create a histogram of the SNRs for a given set of results
Histograms are in log10

Inputs:
    results = list of result dictionaries
    D = telescope diameter (m)
    t = exposure time (s)
    eta = throughput

Outputs:
    List of SNRs for each kernel and wavelength
"""
def calc_SNR_hist(results,D,t,eta):
    SNR_arr = []
    for item in results:
        temp_snr_arr = grab_SNR_per_kernel(item,D,t,eta)
        SNR_arr.append(total_SNR(temp_snr_arr))

    print(np.sum(np.array(SNR_arr)>5))

    plt.figure(1)
    plt.hist(np.log10(SNR_arr),20,histtype="step")
    print(len(SNR_arr))
    return


"""
Create a pandas dataframe of planet SNR as a function of architecture (for a given mode and wavelength)

Inputs:
    mode = mode index [1,2]
    wave = wavelength index [0,1,2]
    D = telescope diameter (m)
    t = exposure time (s)
    eta = throughput
    baseline_lim = (Boolean) Impose baseline constraints?
    per_telescope = (Boolean) SNR per telescope
    extra_data = (Boolean) attach extra data to df

Output:
    Dataframe of SNR as a function of architecture
"""
def create_dataframe(mode,wave_index,D,t,eta,baseline_lim,per_telescope,extra_data):
    prefix = "out_data/avatar_test"
    mode_names = ["Search","Characterisation"]
    arch = [1,3,4,7,8]
    wavelengths = [10,15,18]
    arch_names = ["Bracewell","Three_telescopes","Four_telescopes","Five_telescopes_K1","Five_telescopes_K2"]
    arch_num_tel = [4,3,4,5,5]

    #Baseline limits
    baseline_min = 10
    baseline_max = 600

    #Name of dataframe (for saving?)
    df_name = prefix+"_dataframe_"+mode_names[mode-1]+"_"+str(wavelengths[wave_index])+".df"

    ls_ser = []
    for j,a in enumerate(arch):
        #Create list of dictionaries and sort them
        results = load_results(prefix,a,mode,wavelengths[wave_index])
        results.sort(key=lambda item: item.get("planet_name"))

        SNR_arr = []
        SNR_key = []
        for item in results:
            SNR_key.append(item['planet_name'])

            #Get SNR
            temp_snr_arr = grab_SNR_per_kernel(item,D,t,eta)
            total = total_SNR(temp_snr_arr)

            #Make SNR 0 if outside the allowable baselines
            if baseline_lim:
                baseline = item['baseline (m)']
                if baseline < baseline_min or baseline > baseline_max:
                    total = 0

            #Divide SNR by number of telescopes
            if per_telescope:
                total /= arch_num_tel[j]
            SNR_arr.append(total)
        ls_ser.append(pd.Series(SNR_arr, index=SNR_key, copy=False, name=arch_names[j]))

    #Add spectral type, planet angular separation and distance to the dataframe
    if extra_data:
        type = []
        planet_sep = []
        dist = []
        for item in results:
            type.append(item["star_type"])
            planet_sep.append(item['planet_angle (mas)'])
            dist.append(item["star_distance (pc)"])
        ls_ser.append(pd.Series(type, index=SNR_key, copy=False, name="star_type"))
        ls_ser.append(pd.Series(planet_sep, index=SNR_key, copy=False, name="planet_sep"))
        ls_ser.append(pd.Series(dist, index=SNR_key, copy=False, name="dist"))

    df = pd.concat(ls_ser,axis=1)
    return df


"""
Make a pie chart of the architecture with best SNR over all planets
Inputs:
    mode = mode index [1,2]
    wave = wavelength [10,15,18]
    baseline_lim = (Boolean) Impose baseline constraints?
    per_telescope = (Boolean) SNR per telescope

Output:
    Dataframe of SNR as a function of architecture
"""
def make_SNR_pie_chart(mode,wave_index,baseline_lim,per_telescope):
    D = 2 #2m diameter
    t = 3600 #1hr integration
    eta = 0.1 #10% throughput

    df = create_dataframe(mode,wave_index,D,t,eta,baseline_lim,per_telescope,0)

    maxidx = df.idxmax(axis=1)
    maxidx.value_counts().plot(kind='pie')


def main(arch,mode,wave,D,t,eta):
    prefix = "out_data/avatar_test"
    results = load_results(prefix,arch,mode,wave)
    calc_SNR_hist(results,D,t,eta)
