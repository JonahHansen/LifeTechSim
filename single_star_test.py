import numpy as np
import pickle
from five_nuller_computer import compute
from sim_functions import calc_local_zodiacal_minimum,Spectrograph
from planet_retrieval import RetrievePlanetData as RPD

#Main parameters
sz = 400
fov_scale_factor = 3
mode = 2

first_run = False

spec = Spectrograph(15.091,2.923,15,10)


if first_run:

    planet_path = "PPop/TestPlanetPopulation2.txt"

    star_list = RPD(planet_path,spec)

    f = open("sst_star_list.pkl","wb")

    pickle.dump(star_list, f)

    f.close()

else:

    f = open("sst_star_list.pkl","rb")

    star_list = pickle.load(f)

    f.close()



local_exozodi = calc_local_zodiacal_minimum(spec)

star_data = compute(star_list[0],mode,spec,sz,fov_scale_factor,local_exozodi)
