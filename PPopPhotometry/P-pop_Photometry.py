"""
# =============================================================================
# P-POP PHOTOMETRY
# A photometry tool for P-POP
#
# Authors: Jens Kammerer, Sascha Quanz, Emile Fontanet
# Version: 5.0.0
# Last edited: 26.03.2021
# =============================================================================
#
# P-pop is introduced in Kammerer & Quanz 2018
# (https://ui.adsabs.harvard.edu/abs/2018A%26A...609A...4K/abstract). Please
# cite this paper if you use P-pop for your research.
#
# P-pop makes use of forecaster from Chen & Kipping 2017
# (https://ui.adsabs.harvard.edu/abs/2017ApJ...834...17C/abstract).
"""


# Don't print annoying warnings. Comment out if you want to see them.
import warnings
warnings.filterwarnings('ignore')


# =============================================================================
# IMPORTS
# =============================================================================

# Import your own filters and photometry tools here.
import PhotometryComputer
from Filters import SVO
from Star import Blackbody
from Planet import Thermal, Reflected


# =============================================================================
# SETUP
# =============================================================================

# Select the filters and photometry tools which you want to use here.

# Select the name of the planet population table for which the photometry
# should be computed here.
PathPlanetTable = '../P-pop/TestPlanetPopulation.txt' # str

# Select the filters for which the photometry should be computed here. You can
# simply use the filter names from the Spanish Virtual Observatory
# (http://svo2.cab.inta-csic.es/theory/fps/).
# list of str
SVOids = ['JWST/MIRI.F560W',\
          'JWST/MIRI.F1000W',\
          'JWST/MIRI.F1500W'] # used for LIFE
#SVOids = ['Paranal/SPHERE.ZIMPOL_V',\
#          'Paranal/SPHERE.IRDIS_B_J',\
#          'Paranal/SPHERE.IRDIS_B_H'] # used for HabEx/LUVOIR

# Select the photometry tools to compute the fluxes from the stars and the
# planets as well as their unit and the wavelength range in which the mission
# is operating here.
Sstar = [Blackbody] # list of Star
Splanet = [Thermal, Reflected] # list of Planet
Unit = 'uJy' # micro-Jansky
#Unit = 'ph' # photons per second per square meter
Mission = 'MIR' # use AgeomMIR for reflected light (used for LIFE)
#Mission = 'VIS' # use AgeomVIS for reflected light (used for HabEx/LUVOIR)

# Select whether you want to display summary plots after loading the filters
# and models selected above.
SummaryPlots = True
#SummaryPlots = False
#FigDir = None # if you don't want to save the summary plots
FigDir = 'Figures/' # should end with a slash ("/")
#block = True
block = False


# =============================================================================
# P-POP PHOTOMETRY
# =============================================================================

# Don't modify the following code.
Filters = []
for i in range(len(SVOids)):
    Filters += [SVO.getFilter(SVOids[i], SummaryPlots, FigDir, block)]

PhotComp = PhotometryComputer.PhotometryComputer(PathPlanetTable,
                                                 Filters,
                                                 Sstar,
                                                 Splanet,
                                                 Unit,
                                                 Mission,
                                                 SummaryPlots,
                                                 FigDir,
                                                 block)
PhotComp.Run()
