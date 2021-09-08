import numpy as np
from planet_retrieval import RetrievePlanetData as RPD
from sim_computer import compute

import PPopPhotometry.PhotometryComputer
from PPopPhotometry.Filters import SVO
from PPopPhotometry.Star import Blackbody
from PPopPhotometry.Planet import Thermal, Reflected

#Main parameters
sz = 400
telescope_diameter =
exp_time =
fov_scale_factor = 3 #Need better name - number of phase cycles to
eta =
mode =

#Important derived parameters
telescope_area = np.pi*(telescope_diameter/2)**2

planet_path =

#############################################################################
#Generate photometry

# Select the filters for which the photometry should be computed here. You can
# simply use the filter names from the Spanish Virtual Observatory
# (http://svo2.cab.inta-csic.es/theory/fps/).
# str
SVOid = 'JWST/MIRI.F560W'
#       'JWST/MIRI.F1000W'
#       'JWST/MIRI.F1500W'

# Select the photometry tools to compute the fluxes from the stars and the
# planets as well as their unit and the wavelength range in which the mission
# is operating here.
Sstar = [Blackbody] # list of Star
Splanet = [Thermal, Reflected] # list of Planet
Unit = 'ph' # photons per second per square meter
Mission = 'MIR' # use AgeomMIR for reflected light (used for LIFE)

# Select whether you want to display summary plots after loading the filters
# and models selected above.
SummaryPlots = False
FigDir = None # if you don't want to save the summary plots
block = True

# Don't modify the following code.
filter = SVO.getFilter(SVOid, SummaryPlots, FigDir, block)

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

for star in star_list:
    compute(star,mode,filter,sz,fov_scale_factor,telescope_area,exp_time,eta)
