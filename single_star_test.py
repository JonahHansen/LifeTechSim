import numpy as np
import pickle
from engine.main_computer import compute
from engine.sim_functions import calc_local_zodiacal_minimum,Spectrograph
from engine.planet_retrieval import RetrievePlanetData as RPD

#Main parameters
mode = 2 #1 is search mode, 2 is characterisation mode
spec = Spectrograph(3,18,15,10) #min,max,baseline_wavelength,num_channels
planet_path = "PPop/TestPlanetPopulation2.txt" #Input planet data

if mode == 1:
    sz = 1500
    fov_scale_factor = 5
elif mode == 2:
    sz = 600
    fov_scale_factor = 2


#Run planet retrieval (to convert from PPop file), or load pickled file?
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

#Get local zodi minimum
local_exozodi = calc_local_zodiacal_minimum(spec)

#Choose which architecture and baseline scale factor
from engine.nullers.bracewell import get_nuller_response
base_scale_factor = 0.59

#Which star?
star_index = 2

#run
star_data = compute(star_list[star_index],mode,get_nuller_response,spec,sz,base_scale_factor,fov_scale_factor,local_exozodi)

print(star_data)
