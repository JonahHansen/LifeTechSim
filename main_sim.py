import numpy as np
import pickle
from engine.sim_functions import calc_local_zodiacal_minimum,Spectrograph
from engine.planet_retrieval import RetrievePlanetData as RPD
from engine.main_computer import compute
from itertools import chain
from multiprocessing import Pool
import json
import sys

if len(sys.argv) != 6:
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
"""
architecture = int(sys.argv[1])
mode = int(sys.argv[2])
base_wave = float(sys.argv[3]) #microns
out_file = str(sys.argv[4])
first_run = bool(sys.argv[5])

#####################################################
#Secondary parameters

min_wave = 3 #microns
max_wave = 18 #microns
num_channels = 10

planet_path = "PPop/LifeTechSimTestPlanetPopulation.txt"

number_processes = 28
#####################################################

spec = Spectrograph(min_wave,max_wave,base_wave,num_channels)

#Set architecture
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

else:
    raise Exception("Architecture not recognised")

#Set modes
if mode == 1:
    sz = 1500
    mode_verbose = "Search"
    fov_scale_factor = 5
elif mode == 2:
    sz = 600
    mode_verbose = "Characterisation"
    fov_scale_factor = 2
else:
    raise Exception("Mode not recognised")

#Get star list
if first_run:

    star_list = RPD(planet_path,spec)

    f = open("sst_star_list.pkl","wb")
    pickle.dump(star_list, f)
    f.close()

else:
    f = open("sst_star_list.pkl","rb")
    star_list = pickle.load(f)
    f.close()

#Get local zodi minimum
local_exozodi = calc_local_zodiacal_minimum(spec)

#RUN!!

def worker_func(star):
    return compute(star,mode,get_nuller_response,spec,sz,base_scale_factor,fov_scale_factor,local_exozodi)

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
