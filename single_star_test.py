import numpy as np
import pickle
from five_nuller_computer import compute
from sim_functions import calc_local_zodiacal_minimum
from PPopPhotometry.Filters import SVO

#Main parameters
sz = 400
fov_scale_factor = 3 #Need better name - number of phase cycles to
mode = 2

first_run = False

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
filters = [filter]

if first_run:

    from planet_retrieval import RetrievePlanetData as RPD
    import PPopPhotometry.PhotometryComputer as PhotometryComputer
    from PPopPhotometry.Star import Blackbody
    from PPopPhotometry.Planet import Thermal, Reflected

    planet_path = "PPop/TestPlanetPopulation2.txt"


    #############################################################################
    #Generate photometry

    # Select the photometry tools to compute the fluxes from the stars and the
    # planets as well as their unit and the wavelength range in which the mission
    # is operating here.
    Sstar = [Blackbody] # list of Star
    Splanet = [Thermal, Reflected] # list of Planet
    Unit = 'ph' # photons per second per square meter
    Mission = 'MIR' # use AgeomMIR for reflected light (used for LIFE)

    PhotComp = PhotometryComputer.PhotometryComputer(planet_path,
                                                     filters,
                                                     Sstar,
                                                     Splanet,
                                                     Unit,
                                                     Mission,
                                                     SummaryPlots,
                                                     FigDir,
                                                     block)
    PhotComp.Run()

    phot_path = planet_path.split(".")[0]+'_'+SVOid.split('/')[1]+'.txt'

    ########################################################################

    star_list = RPD(planet_path,phot_path)

    f = open("sst_star_list.pkl","wb")

    pickle.dump(star_list, f)

    f.close()

else:

    f = open("sst_star_list.pkl","rb")

    star_list = pickle.load(f)

    f.close()

local_exozodi = calc_local_zodiacal_minimum(filter)

star_data = compute(star_list[0],mode,filter,sz,fov_scale_factor,local_exozodi)
