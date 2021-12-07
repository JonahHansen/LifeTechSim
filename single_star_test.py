import numpy as np
import pickle
from engine.main_computer import compute
from engine.sim_functions import calc_local_zodiacal_minimum,Spectrograph
from engine.planet_retrieval import RetrievePlanetData as RPD
from snr_calculator import total_SNR, grab_SNR_per_kernel

#Main parameters
mode = 2 #1 is search mode, 2 is characterisation mode
architecture = 2
main_wave = 15
spec = Spectrograph(3,18,main_wave,10) #min,max,baseline_wavelength,num_channels
planet_path = "PPop/TestPlanetPopulation2.txt" #Input planet data

if mode == 1:
    sz = 1500
    fov_scale_factor = 5
elif mode == 2:
    sz = 800
    fov_scale_factor = main_wave/(spec.channel_centres[0]*1e6) + 0.1


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

#Set architecture, and define the baseline scale factor
if architecture == 1:
    from engine.nullers.bracewell import get_nuller_response
    architecture_verbose = "Bracewell four telescope nuller"
    base_scale_factor = 0.590

elif architecture == 2:
    from engine.nullers.linear import get_nuller_response
    architecture_verbose = "Linear four telescope nuller"
    base_scale_factor = 0.6

elif architecture == 3:
    from engine.nullers.three_telescopes import get_nuller_response
    architecture_verbose = "Three telescope kernel nuller"
    base_scale_factor = 0.666

elif architecture == 4:
    from engine.nullers.four_telescopes import get_nuller_response
    architecture_verbose = "Four telescope kernel nuller, optimised for K1"
    base_scale_factor = 0.4

elif architecture == 5:
    from engine.nullers.linear_assymmetric import get_nuller_response
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


#Which star?
star_index = 2

#run
star_data = compute(star_list[star_index],mode,get_nuller_response,spec,sz,base_scale_factor,fov_scale_factor,local_exozodi)

snr_arr=[]
for dict in star_data:
    if architecture == (4 or 6):
        kernel_snr = grab_SNR_per_kernel(dict,2,3600,0.1,0.5)
    else:
        kernel_snr = grab_SNR_per_kernel(dict,2,3600,0.1,1)
    snr_arr.append(total_SNR(kernel_snr))

print(snr_arr)
