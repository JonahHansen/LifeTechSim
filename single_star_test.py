import numpy as np
import pickle
from engine.five_nuller_computer import compute
from engine.sim_functions import calc_local_zodiacal_minimum,Spectrograph
from planet_retrieval import RetrievePlanetData as RPD

#Main parameters
mode = 1 #1 is search mode, 2 is characterisation mode
spec = Spectrograph(3,18,15,10) #min,max,baseline_wavelength,num_channels
planet_path = "PPop/TestPlanetPopulation2.txt"

if mode == 1:
    sz = 1500
    fov_scale_factor = 5
elif mode == 2:
    sz = 600
    fov_scale_factor = 2

first_run = False

if first_run:

    star_list = RPD(planet_path,spec)

    f = open("sst_star_list.pkl","wb")
    pickle.dump(star_list, f)
    f.close()

else:
    f = open("sst_star_list.pkl","rb")
    star_list = pickle.load(f)
    f.close()

local_exozodi = calc_local_zodiacal_minimum(spec)

star_data = compute(star_list[5],mode,spec,sz,fov_scale_factor,local_exozodi)
