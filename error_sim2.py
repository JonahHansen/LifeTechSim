import numpy as np
import pickle
from engine.sim_functions import calc_local_zodiacal_minimum,Spectrograph
from engine.planet_retrieval import Planet,Star
from engine.main_computer import compute
from itertools import chain
from multiprocessing import Pool
import json
import sys

dR_scale = float(sys.argv[1]) #Reflectance error in beamsplitter
out_file = str(sys.argv[2])

base_wave = 18
#####################################################
#Secondary parameters

min_wave = 4 #microns
max_wave = 19 #microns
num_channels = 50

#####################################################

#Set up the spectral parameters
master_spec = Spectrograph(min_wave,max_wave,base_wave,num_channels)

central_waves = min_wave/2*np.exp(1/3*(np.log(max_wave)-np.log(min_wave))*np.array([1,2,3]))*(1+np.exp(-1/3*(np.log(max_wave)-np.log(min_wave))))

new_specs = []
for i in master_spec.channel_borders:
    new_spec = Spectrograph(i,i+master_spec.dlambda,base_wave,1)
    new_specs.append(new_spec)

main_dict_ls = []

from engine.nullers.five_telescopes_err_2 import get_nuller_response
architecture_verbose = "Five telescope kernel nuller, optimised for diagonal telescopes (K2 alt)"
base_scale_factor = 1.028 #= approx 1.03*0.619 (where 0.619 is the conversion between a side and diagonal of a pentagon)

sz = 1500
mode_verbose = "Search"
fov_scale_factor = 5

dphi_scale = dR_scale
number_processes = 1 #parallelise?

    ###########################################################################
for spec in new_specs:
    #Tau boo coordinates (ecliptic latitude of 26 degress)
    ra = 206.81560195
    dec = 17.45690446
    z = 1 #exozodi (same as solar system)

    #Generate stars of given types
    def Beta_Pic(dist):
        name = "Beta Pic analogue"
        stype = "A6V"
        rad = 1.8
        teff = 8052
        mass = 1.75
        return Star(name,1,dist,stype,rad,teff,mass,ra,dec,spec,z)

    def Tau_Boo(dist):
        name = "Tau Boo analogue"
        stype = "F7V"
        rad = 1.42
        teff = 6399
        mass = 1.39
        return Star(name,2,dist,stype,rad,teff,mass,ra,dec,spec,z)

    def Sun(dist):
        name = "Solar analogue"
        stype = "G2V"
        rad = 1
        teff = 5772
        mass = 1
        return Star(name,3,dist,stype,rad,teff,mass,ra,dec,spec,z)

    def Eps_Eri(dist):
        name = "Epsilon Eri analogue"
        stype = "K2V"
        rad = 0.735
        teff = 5084
        mass = 0.82
        return Star(name,4,dist,stype,rad,teff,mass,ra,dec,spec,z)

    def Prox_Cen(dist):
        name = "Proxima Cen analogue"
        stype = "M5V"
        rad = 0.15
        teff = 3042
        mass = 0.12
        return Star(name,5,dist,stype,rad,teff,mass,ra,dec,spec,z)

    #Helper function to generate Earth-twin planets
    def myPlanet(star,num,a):

        #Earth twin
        PRad = 1
        PMass = 1
        Temp = 300

        Ageom = 0.1 #Rough estimate?
        AngSep = a/star.Dist
        lam_ref = 0.318 #From PPop, assuming face on orbit (inc = 0)

        return Planet(star,0,star.SNumber,num,PRad,PMass,365,0,0,0,0,0,0,0,Ageom,a,a,AngSep,0,0,lam_ref,Temp,spec)

    #Give each star one planet in the middle of the HZ
    def append_planet_list(star):
        star.Planets = [myPlanet(star,2,star.HZMid)]
        return star

    dists = np.linspace(1,20,20)

    #Make the list of stars at given distances
    star_list = []
    for d in dists:
        star_list.append(append_planet_list(Sun(d)))
        star_list.append(append_planet_list(Eps_Eri(d)))
        star_list.append(append_planet_list(Prox_Cen(d)))


    #Make errors
    dphi = np.zeros(10)
    dR = np.zeros(10)

    phases = [-2.055,-2.840,0.810,2.997,-0.494,-1.571]
    phase_chops = np.abs(np.array([2*phases,2*np.pi-2*np.abs(phases)]))
    min_phase_chops = phase_chops.T[np.arange(len(phase_chops.T)),np.argmin(phase_chops,axis=0)]

    wave = spec.mean/1.6
    centre = np.min(np.abs(channel_centres-wave))+wave
    phase_chop_errs = np.abs(min_phase_chops*wave/centre - min_phase_chops)


    dphi[2] = phase_chop_errs[0]+np.sign(np.random.random()*2-1)*dphi_scale
    dphi[4] = phase_chop_errs[1]+np.sign(np.random.random()*2-1)*dphi_scale
    dphi[5] = phase_chop_errs[2]+np.sign(np.random.random()*2-1)*dphi_scale
    dphi[7] = phase_chop_errs[3]+np.sign(np.random.random()*2-1)*dphi_scale
    dphi[8] = phase_chop_errs[4]+np.sign(np.random.random()*2-1)*dphi_scale
    dphi[9] = phase_chop_errs[5]+np.sign(np.random.random()*2-1)*dphi_scale

    dR[2] = np.sign(np.random.random()*2-1)*dR_scale
    dR[4] = np.sign(np.random.random()*2-1)*dR_scale
    dR[5] = np.sign(np.random.random()*2-1)*dR_scale
    dR[7] = np.sign(np.random.random()*2-1)*dR_scale
    dR[8] = np.sign(np.random.random()*2-1)*dR_scale
    dR[9] = np.sign(np.random.random()*2-1)*dR_scale

    def response_func(baseline,fov,sz,base_wavelength):
        return get_nuller_response(dphi,dR,baseline,fov,sz,base_wavelength)

    ###########################################################################

    #Get local zodi minimum
    local_exozodi = calc_local_zodiacal_minimum(spec)

    #RUN!!

    #Multiprocessing
    def worker_func(star):
        return compute(star,1,response_func,spec,sz,base_scale_factor,fov_scale_factor,local_exozodi)

    pool = Pool(processes=number_processes)
    ls_star_data = pool.map(worker_func,star_list)
    pool.close()

    #Make into one list of dictionaries
    dict_ls = list(chain.from_iterable(ls_star_data))

    #Funtion to round sig figs (for numerical readibility)
    def round_sig_figs(x, p):
        x_positive = np.where(np.isfinite(x) & (x != 0), np.abs(x), 10**(p-1))
        mags = 10 ** (p - 1 - np.floor(np.log10(x_positive)))
        return np.round(x * mags) / mags

    #Header data, plus results
    main_dict = {
        "Architecture":architecture_verbose,
        "Mode":mode_verbose,
        "baseline_wavelength (microns)":round_sig_figs(spec.baseline_wave*1e6,5),
        "sz":sz,
        "fov_scale_factor":fov_scale_factor,
        "min_wavelength (microns)":round_sig_figs(spec.wave_min*1e6,5),
        "max_wavelength (microns)":round_sig_figs(spec.wave_max*1e6,5),
        "num_channels":spec.n_channels,
        "channel_widths (microns)":round_sig_figs(spec.dlambda*1e6,5),
        "channel_centres (microns)":round_sig_figs(spec.channel_centres*1e6,5),
        "results":dict_ls
    }

    main_dict_ls.append(main_dict)

main_dict = {
    "Architecture":architecture_verbose,
    "Mode":mode_verbose,
    "baseline_wavelength (microns)":round_sig_figs(master_spec.baseline_wave*1e6,5),
    "sz":sz,
    "fov_scale_factor":fov_scale_factor,
    "min_wavelength (microns)":round_sig_figs(master_spec.wave_min*1e6,5),
    "max_wavelength (microns)":round_sig_figs(master_spec.wave_max*1e6,5),
    "num_channels":master_spec.n_channels,
    "channel_widths (microns)":round_sig_figs(master_spec.dlambda*1e6,5),
    "channel_centres (microns)":round_sig_figs(master_spec.channel_centres*1e6,5),
    "results":main_dict_ls
}

#Needed for writing JSON
class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)

#Write file
fout = open(out_file,"w")
json.dump(main_dict,fout,cls=NpEncoder,indent=2)
fout.close()
