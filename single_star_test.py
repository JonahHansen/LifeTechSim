import numpy as np
import pickle
from five_nuller_computer import compute
from sim_functions import calc_local_zodiacal_minimum
from PPopPhotometry.Filters import SVO

#Main parameters
sz = 400
fov_scale_factor = 3 #Need better name - number of phase cycles to
mode =

first_run = True

# Select the filters for which the photometry should be computed here. You can
# simply use the filter names from the Spanish Virtual Observatory
# (http://svo2.cab.inta-csic.es/theory/fps/).
# str
SVOid = 'JWST/MIRI.F1500W'
#       'JWST/MIRI.F1000W'

# Select whether you want to display summary plots after loading the filters
# and models selected above.
SummaryPlots = False
FigDir = None # if you don't want to save the summary plots
block = True

# Don't modify the following code.
filter = SVO.getFilter(SVOid, SummaryPlots, FigDir, block)

if first_run:

    from planet_retrieval import RetrievePlanetData as RPD
    import PPopPhotometry.PhotometryComputer
    from PPopPhotometry.Star import Blackbody
    from PPopPhotometry.Planet import Thermal, Reflected

    planet_path =


    #############################################################################
    #Generate photometry

    # Select the photometry tools to compute the fluxes from the stars and the
    # planets as well as their unit and the wavelength range in which the mission
    # is operating here.
    Sstar = [Blackbody] # list of Star
    Splanet = [Thermal, Reflected] # list of Planet
    Unit = 'ph' # photons per second per square meter
    Mission = 'MIR' # use AgeomMIR for reflected light (used for LIFE)

    PhotComp = PhotometryComputer.PhotometryComputer(PathPlanetTable,
                                                     filter,
                                                     Sstar,
                                                     Splanet,
                                                     Unit,
                                                     Mission,
                                                     SummaryPlots,
                                                     FigDir,
                                                     block)
    PhotComp.Run()

    phot_path = planet_path+'_'+SVOid.split('/')[1]+'.txt'

    ########################################################################

    star_list = RPD(planet_path,phot_path)

    pickle.dump(star_list, "sst_star_list.pkl")

else:

    star_list = pickle.load("sst_star_list.pkl")

local_exozodi = calc_local_zodiacal_minimum(filter)

star_data = compute(star_list[0],mode,filter,sz,fov_scale_factor,local_exozodi)