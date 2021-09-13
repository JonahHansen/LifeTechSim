import numpy as np
from itertools import chain
from planet_retrieval import RetrievePlanetData as RPD
from sim_computer import compute
from astropy.table import Table

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
output_path =

# Select the filters for which the photometry should be computed here. You can
# simply use the filter names from the Spanish Virtual Observatory
# (http://svo2.cab.inta-csic.es/theory/fps/).
# str
SVOid = 'JWST/MIRI.F1500W'
#       'JWST/MIRI.F1000W'

#############################################################################
#Generate photometry

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

ls_star_data = []
for star in star_list:
    star_data = compute(star,mode,filter,sz,fov_scale_factor,telescope_area,exp_time,eta)
    ls_star_data.append(star_data)

dict_ls = list(chain.from_iterable(ls_star_data))
Fits_table = Table(rows=dict_ls)

Fits_Table.meta["Size"] = sz
Fits_Table.meta["Architecture"] = ""
Fits_Table.meta["Diameter"] = telescope_diameter
Fits_Table.meta["Scan_Mode"] = ""
Fits_Table.meta["Filter"] = filter.Name
Fits_Table.meta["FOV_scale"] = fov_scale_factor

Fits_table.write(output_path,format='fits')
