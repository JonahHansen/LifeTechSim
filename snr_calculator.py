import json
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
plt.ion()


def SNR(signal,shot,leakage,zodiacal,exozodiacal):
    sig = signal
    noise = 2*shot + 2*leakage + 2*zodiacal + 2*exozodiacal
    return sig/np.sqrt(noise)

def SNR_collapse(SNR_array):
    return np.sqrt(np.sum(SNR_array**2))

def grab_SNR_per_kernel(dict,D,t,eta):

    A = np.pi*D**2/4

    signal = np.array(dict["signal (ph/s/m2)"])*A*t*eta
    shot = np.array(dict["shot (ph/s/m2)"])*A*t*eta
    leakage = np.array(dict["leakage (ph/s/m2)"])*A*t*eta
    exozodiacal = np.array(dict["exozodiacal (ph/s/m2)"])*A*t*eta
    zodiacal = np.array(dict["zodiacal (ph/s)"])*t*eta

    return SNR(signal, shot, leakage, zodiacal, exozodiacal)

def total_SNR(SNR_array):
    return SNR_collapse(SNR_array.flatten())

def cost(N,D,lam,T=20):
    return N*300*D**1.7*lam**(-0.5)*T**(-0.25)

def load_results(prefix,arch,mode,wave):
    filename = prefix+"_"+str(arch)+"_"+str(mode)+"_"+str(wave)+".json"
    fin = open(filename,"r")
    data = json.load(fin)
    fin.close()
    return data["results"]

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

def main(arch,mode,wave,D,t,eta):
    prefix = "out_data/avatar_test"
    results = load_results(prefix,arch,mode,wave)
    calc_SNR_hist(results,D,t,eta)

def create_dataframe(mode,wave_index,D,t,eta,baseline_lim,per_telescope):
    prefix = "out_data/avatar_test"
    mode_names = ["Search","Characterisation"]
    arch = [1,3,4,7,8]
    wavelengths = [10,15,18]
    arch_names = ["Bracewell","Three_telescopes","Four_telescopes","Five_telescopes_K1","Five_telescopes_K2"]
    arch_num_tel = [4,3,4,5,5]

    baseline_min = 10
    baseline_max = 600

    df_name = prefix+"_dataframe_"+mode_names[mode-1]+"_"+str(wavelengths[wave_index])+".df"

    ls_ser = []
    for j,a in enumerate(arch):
        results = load_results(prefix,a,mode,wavelengths[wave_index])
        results.sort(key=lambda item: item.get("planet_name"))
        SNR_arr = []
        SNR_key = []
        for item in results:
            SNR_key.append(item['planet_name'])
            temp_snr_arr = grab_SNR_per_kernel(item,D,t,eta)
            total = total_SNR(temp_snr_arr)
            if baseline_lim:
                baseline = item['baseline (m)']
                if baseline < baseline_min or baseline > baseline_max:
                    total = 0
            if per_telescope:
                total /= arch_num_tel[j]
            SNR_arr.append(total)
        ls_ser.append(pd.Series(SNR_arr, index=SNR_key, copy=False, name=arch_names[j]))

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

#maxidx = df.idxmax(axis=1)
#maxidx.value_counts().plot(kind='pie')
