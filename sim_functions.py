import numpy as np


def limb_darkening(r):
    mu = np.sqrt(1-r**2)
    return 1-0.47*(1-mu)-0.23*(1-mu)**2


def get_SNR(signal,leakage,zodiacal,N,area,exp_time,eta):

    C_per_t = area*exp_time*eta
    #Signal
    signal_per_t = np.sum(signal*C_per_t)
    #Noise contributions: shot, leakage and zodiacal background
    noise_per_t = np.sqrt(np.sum(signal*C_per_t) + np.sum(leakage*C_per_t) + zodiacal*C_per_t)
    SNR_per_t = signal_per_t/noise_per_t


    C_tot = N*area*exp_time*eta
    #Signal
    signal_tot = np.sum(signal*C_tot)
    #Noise contributions: shot, leakage and zodiacal background
    noise_tot = np.sqrt(np.sum(signal*C_tot) + np.sum(leakage*C_tot) + zodiacal*C_tot)
    SNR_tot = signal_tot/noise_tot

    return SNR_per_t, SNR_tot
