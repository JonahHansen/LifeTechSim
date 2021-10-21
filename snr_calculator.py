import json
import numpy as np

filename = "avatar_test.json"


def SNR(signal,leakage,zodiacal,exozodiacal):
    sig = signal
    noise = 2*signal + 2*leakage + zodiacal + exozodiacal
    return sig/np.sqrt(noise)


def grab_SNR_per_telescope(dict,D,t,eta):

    A = np.pi*D**2/4

    signal_k1 = np.array(dict["signal_k1 (ph/s/m2)"])*A*t*eta
    leakage_1 = np.array(dict["leakage_1 (ph/s/m2)"])*A*t*eta
    exozodiacal_1 = np.array(dict["exozodiacal_1 (ph/s/m2)"])*A*t*eta
    signal_k2 = np.array(dict["signal_k2 (ph/s/m2)"])*A*t*eta
    leakage_2 = np.array(dict["leakage_2 (ph/s/m2)"])*A*t*eta
    exozodiacal_2 = np.array(dict["exozodiacal_2 (ph/s/m2)"])*A*t*eta
    zodiacal = np.array(dict["zodiacal (ph/s)"])*t*eta

    snr_1 = SNR(signal_k1,leakage_1,zodiacal,exozodiacal_1)
    snr_2 = SNR(signal_k2,leakage_2,zodiacal,exozodiacal_2)

    return snr_1,snr_2

def grab_SNR_total(dict,N,D,t,eta):
    A = N*np.pi*D**2/4

    signal_k1 = np.array(dict["signal_k1 (ph/s/m2)"])*A*t*eta
    leakage_1 = np.array(dict["leakage_1 (ph/s/m2)"])*A*t*eta
    exozodiacal_1 = np.array(dict["exozodiacal_1 (ph/s/m2)"])*A*t*eta
    signal_k2 = np.array(dict["signal_k2 (ph/s/m2)"])*A*t*eta
    leakage_2 = np.array(dict["leakage_2 (ph/s/m2)"])*A*t*eta
    exozodiacal_2 = np.array(dict["exozodiacal_2 (ph/s/m2)"])*A*t*eta
    zodiacal = np.array(dict["zodiacal (ph/s)"])*t*eta

    snr_1 = SNR(signal_k1,leakage_1,zodiacal,exozodiacal_1)
    snr_2 = SNR(signal_k2,leakage_2,zodiacal,exozodiacal_2)

    return snr_1, snr_2

def SNR_collapse_wavelengths(snr_array):
    return np.sqrt(np.sum(snr_array**2))

fin = open(filename,"r")

data = json.load(fin)

fin.close()

results = data["results"]

x = results[0]
