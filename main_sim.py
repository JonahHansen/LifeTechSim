import numpy as np
from planet_retrieval import RetrievePlanetData as RPD
from sim_computer import compute


#Main parameters
wavelength = 10e-6
sz = 400
telescope_diameter =
exp_time =
fov_scale_factor = 3 #Need better name - number of phase cycles to
eta =
mode =

#Important derived parameters
telescope_area = np.pi*(telescope_diameter/2)**2

planet_path =
phot_path =

star_list = RPD(planet_path,phot_path)

for star in star_list:
    compute(star,mode,wavelength,sz,fov_scale_factor,telescope_area,exp_time,eta)
