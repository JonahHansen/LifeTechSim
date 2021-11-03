import numpy as np
import pickle
from engine.sim_functions import calc_local_zodiacal_minimum,Spectrograph
from engine.planet_retrieval import Planet,Star
from engine.main_computer import compute
from itertools import chain
from multiprocessing import Pool
import json
import sys

if len(sys.argv) != 4:
    raise Exception("Wrong number of arguments")

#####################################################
#Main parameters, passed as commandline arguments

"""
Architecture
 1 = Bracewell
 2 = Linear nuller
 3 = Three way kernel nuller
 4 = Four way kernel nuller (Kernel 1)
 5 = Four way kernel nuller (Kernel 2)
 6 = Four way kernel nuller (Kernel 3)
 7 = Five way kernel nuller (Kernel 1)
 8 = Five way kernel nuller (Kernel 2)

Mode
 1 = Search (optimises habitible zone of star, array spins)
 2 = Charaterisation (optimised planet position, no spinning)

Base_wave = wavelength to base the optimisation on

out_file = output filename

first_run = run planet retrieval (to convert from PPop file) or import from pickle?
"""
architecture = int(sys.argv[1])
base_wave = float(sys.argv[2]) #microns
out_file = str(sys.argv[3])

#####################################################
#Secondary parameters

min_wave = 3 #microns
max_wave = 18 #microns
num_channels = 10

#input planet data
planet_path = "PPop/LifeTechSimTestPlanetPopulation.txt"

number_processes = 28 #parallelise?
#####################################################

#Set up the spectral parameters
spec = Spectrograph(min_wave,max_wave,base_wave,num_channels)

#Set architecture, and define the baseline scale factor
if architecture == 1:
    from engine.nullers.bracewell import get_nuller_response
    architecture_verbose = "Bracewell four telescope nuller"
    base_scale_factor = 0.590

elif architecture == 2:
    from engine.nullers.linear import get_nuller_response
    architecture_verbose = "Linear four telescope nuller"

elif architecture == 3:
    from engine.nullers.three_telescopes import get_nuller_response
    architecture_verbose = "Three telescope kernel nuller"
    base_scale_factor = 0.666

elif architecture == 4:
    from engine.nullers.four_telescopes import get_nuller_response
    architecture_verbose = "Four telescope kernel nuller, optimised for K1"
    base_scale_factor = 0.4

elif architecture == 5:
    from engine.nullers.four_telescopes import get_nuller_response
    architecture_verbose = "Four telescope kernel nuller, optimised for K2"
    base_scale_factor = 0.4

elif architecture == 6:
    from engine.nullers.four_telescopes import get_nuller_response
    architecture_verbose = "Four telescope kernel nuller, optimised for K3"
    base_scale_factor = 0.4

elif architecture == 7:
    from engine.nullers.five_telescopes import get_nuller_response
    architecture_verbose = "Five telescope kernel nuller, optimised for adjacent telescopes (K1)"
    base_scale_factor = 1.028

elif architecture == 8:
    from engine.nullers.five_telescopes import get_nuller_response
    architecture_verbose = "Five telescope kernel nuller, optimised for diagonal telescopes (K2)"
    base_scale_factor = 0.660 #= approx 1.03*0.619 (where 0.619 is the conversion between a side and diagonal of a pentagon)

elif architecture == 9: #THESE ARE ALTERNATIVES - MAY BE BETTER!
    from engine.nullers.five_telescopes import get_nuller_response
    architecture_verbose = "Five telescope kernel nuller, optimised for adjacent telescopes (K1 alt)"
    base_scale_factor = 2.67

elif architecture == 10:
    from engine.nullers.five_telescopes import get_nuller_response
    architecture_verbose = "Five telescope kernel nuller, optimised for diagonal telescopes (K2 alt)"
    base_scale_factor = 1.68 #= approx 1.03*0.619 (where 0.619 is the conversion between a side and diagonal of a pentagon)


else:
    raise Exception("Architecture not recognised")

sz = 1500
mode_verbose = "Search"
fov_scale_factor = 5

###########################################################################

#Tau boo coordinates (ecliptic latitude of 26 degress)
ra = 206.81560195
dec = 17.45690446
z = 1 #exozodi (same as solar system)

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
    return Star(name,1,dist,stype,rad,teff,mass,ra,dec,spec,z)

def Sun(dist):
    name = "Solar analogue"
    stype = "G2V"
    rad = 1
    teff = 5772
    mass = 1
    return Star(name,1,dist,stype,rad,teff,mass,ra,dec,spec,z)

def Eps_Eri(dist):
    name = "Epsilon Eri analogue"
    stype = "K2V"
    rad = 0.735
    teff = 5084
    mass = 0.82
    return Star(name,1,dist,stype,rad,teff,mass,ra,dec,spec,z)

def Prox_Cen(dist):
    name = "Proxima Cen analogue"
    stype = "M5V"
    rad = 0.15
    teff = 3042
    mass = 0.12
    return Star(name,1,dist,stype,rad,teff,mass,ra,dec,spec,z)

def myPlanet(star,num,a):

    #Earth twin
    PRad = 1
    PMass = 1
    Temp = 300

    Ageom = 0.1 #Rough estimate?
    AngSep = a/star.Dist
    lam_ref = 0.318 #From PPop, assuming face on orbit (inc = 0)

    return Planet(star,0,star.SNumber,num,PRad,PMass,365,0,0,0,0,0,0,0,Ageom,a,a,AngSep,0,0,lam_ref,Temp,spec)

def append_planet_list(star):
    planet_inner = myPlanet(star,1,star.HZIn)
    planet_mid = myPlanet(star,2,star.HZMid)
    planet_outer = myPlanet(star,3,star.HZOut)

    star.Planets = [planet_inner,planet_mid,planet_outer]
    return star


dists = np.linspace(1,20,100)

star_list = []
for d in dists:
    star_list.append(append_planet_list(Beta_Pic(d)))
    star_list.append(append_planet_list(Tau_Boo(d)))
    star_list.append(append_planet_list(Sun(d)))
    star_list.append(append_planet_list(Eps_Eri(d)))
    star_list.append(append_planet_list(Prox_Cen(d)))


###########################################################################

#Get local zodi minimum
local_exozodi = calc_local_zodiacal_minimum(spec)

#RUN!!

#Multiprocessing
def worker_func(star):
    return compute(star,1,get_nuller_response,spec,sz,base_scale_factor,fov_scale_factor,local_exozodi)

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
    "planet_population":planet_path,
    "sz":sz,
    "fov_scale_factor":fov_scale_factor,
    "min_wavelength (microns)":round_sig_figs(spec.wave_min*1e6,5),
    "max_wavelength (microns)":round_sig_figs(spec.wave_max*1e6,5),
    "num_channels":spec.n_channels,
    "channel_widths (microns)":round_sig_figs(spec.dlambda*1e6,5),
    "channel_centres (microns)":round_sig_figs(spec.channel_centres*1e6,5),
    "results":dict_ls
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
