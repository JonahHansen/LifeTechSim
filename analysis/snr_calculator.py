import json
import numpy as np
import matplotlib.pyplot as plt
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
    zod_fac = factor to multiply zodiacal light by (for four telescope kernel nuller). Generally either 1 or 0.5
    const_total_area = whether to keep total area of array the same across different numbers
        of telescopes. Assumes that the telescope diameter is attributed to a four telescope array
    num_telescopes = number of telescopes in array (only used in const_total_area)
    exozod_fac = factor to multiply exozodiacal light by
    stellar_fac = factor to multiply stellar leakage light by

Outputs:
    List of SNRs for each kernel and wavelength
"""
def grab_SNR_per_kernel(dict,D,t,eta,zod_fac=1,const_total_area=False,num_telescopes=4,exozod_fac=1,stellar_fac=1):

    #Conserve total area
    if const_total_area:
        D *= np.sqrt(4/num_telescopes)
    A = np.pi*D**2/4

    signal = np.array(dict["signal"])*A*t*eta
    shot = np.array(dict["shot"])*A*t*eta
    leakage = np.array(dict["leakage"])*A*t*eta*stellar_fac
    exozodiacal = np.array(dict["exozodiacal"])*A*t*eta*exozod_fac
    zodiacal = np.array(dict["zodiacal"])*t*eta*zod_fac

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
Helper function to get the total summed SNR for a given set of results over all kernels and wavelengths

Inputs:
    dict = input result dictionary of a planet (from the main simulation)
    D = telescope diameter (m)
    t = exposure time (s)
    eta = throughput
    zod_fac = factor to multiply zodiacal light by (for four telescope kernel nuller). Generally either 1 or 0.5
    const_total_area = whether to keep total area of array the same across different numbers
        of telescopes. Assumes that the telescope diameter is attributed to a four telescope array
    num_telescopes = number of telescopes in array (only used in const_total_area)

Output:
    Total SNR
"""
def total_SNR_from_dict(dict,D,t,eta,zod_fac=1,const_total_area=False,num_telescopes=4):
    snr_arr = grab_SNR_per_kernel(dict,D,t,eta,zod_fac)
    return total_SNR(snr_arr)


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
def calc_SNR_hist(results,D,t,eta,zod_fac=1,const_total_area=False,num_telescopes=4):
    SNR_arr = []
    for item in results:
        SNR_arr.append(total_SNR_from_dict(item,D,t,eta,zod_fac,const_total_area,num_telescopes))

    print(np.sum(np.array(SNR_arr)>7))

    plt.figure(1)
    plt.hist(np.log10(SNR_arr),20,histtype="step")
    print(len(SNR_arr))
    return
