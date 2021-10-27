import json
import numpy as np

filename = "avatar_test.json"


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
    return SNR_collapse(np.flatten(SNR_array))

def cost(N,D,lambda,T=20):
    return N*300*D**1.7*lambda**(-0.5)*T**(-0.25)

fin = open(filename,"r")

data = json.load(fin)

fin.close()

results = data["results"]

x = results[0]
