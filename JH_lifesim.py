"""
Some code for dealing with the 5 aperture nuller.

For leakage calculations, we need the mean square stellar radius and mean 4th power of
stellar radius for a standard limb-darkening law:

I = I_0 (1 - 0.47(1-mu) - 0.23(1-mu)**2)



Zodiacal light has an optical depth of 1e-7 spread over the sky.
A 290K background has a flux of 2.8e14 Jy/Sr,
2*6.626e-34*nu**3/3e8**2/(np.exp(6.626e-34*nu/1.38e-23/290) - 1)*1e26
... which is 1e7 times higher than the COBE value.
OR
25 MJy/sr
Vega is 40 Jy
Gives 12 magnitudes per square arcsec.

tau Ceti is mag 1.6

For a single nulling baseline with the planet at maximum, the SNR is:
p_tel*2/sqrt(b_tel + stellar_leakage)
b_tel is background per telescope


For a Bracewell nuller, transmission to a single output, with
\alpha_p = lambda/2B
with alpha_p being half fringe spacing in radians

Transmission function of nuller for nulled output
T = 1-cos(2 pi \alpha B/lambda)
  = 1-cos(pi \alpha / \alpha_p)
  = 0.5 \pi^2 (\alpha / \alpha_p)^2 (approxmiation for second order nuller)

The integral of I times T in terms of \alpha/\alpha_p=x is:

\int T(x) I(r)^2 r dr / \int I(r)^2 r dr
\int T(r*cos(theta)) I(r)^2 r dr dtheta / 2*pi \int I(r)^2 r dr

stellar_leakage = integral * s_tel
s_tel is star flux per telescope


"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from astropy import constants as const
import sys
plt.ion()

if not '..' in sys.path:
    sys.path.insert(0,'..')
from opticstools import knull
from opticstools import azimuthalAverage

#Constants in MKS units
sigma_sb = const.sigma_sb.value #Stefan-Boltzmann constant
h = const.h.value #Planck constant
c = const.c.value #Speed of light
k_B = const.k_B.value #Boltzmann constant
R_sol = const.R_sun.value #Solar radius in m
L_sol = const.L_sun.value #Solar luminosity in W
rad2mas = np.degrees(1)*3600e3 #Number of milliarcsec in one radian
au = const.au.value #AU in m
pc = const.pc.value #parsec in m

#Main parameters
wavelength = 10e-6
d_wavelength = 1e-6
sz = 400
number_telescopes =
telescope_diameter =

#Important derived parameters
telescope_area = np.pi*(telescope_diameter/2)**2
total_area = number_telescopes*telescope_area


"""
Return a function for the spectral flux density as a function of wavelength
INPUTS:
    T = Temperature of star in K
    R = Radius of star in R_sol
    d = Distance to star in pc
OUTPUTS:
    Planck function as a function of wavelength (photons/m^2/s/m)
"""
def Planck_function(T,R,d):
    R = R*R_sol #From solar radii to m
    d = d*pc #from parsec to meters
    def func(lam): #lam in m
        B = 2*np.pi*c/(lam**4)/(np.exp(h*c/(lam*k_B*T))-1)*R**2/d**2 #photons/m^2/s/m
        return B
    return func


"""
Return flux density of a star in a given bandpass from Planck's law
INPUTS:
    T = Temperature of star in K
    R = Radius of star in R_sol
    d = Distance to star in pc
    wave = Central wavelength in m
    d_wave = Filter bandpass in m
OUTPUTS:
    Flux density of star in photons/m^2/s
"""
def flux_density(T,R,d,wave,d_wave):
    planck_func = Planck_function(T,R,d)
    I = quad(planck_func,wave-d_wave/2.,wave+d_wave/2.)
    return I[0]


"""
Return luminosity in solar luminosity (from Stefan-Boltzmann law)
INPUTS:
    T = Temperature of star in K
    R = Radius of star in R_sol
OUTPUTS:
    Luminosity of star in L_sol
"""
def Luminosity(T,R):
    Lum = 4*np.pi*(R*R_sol)**2*sigma_sb*T**4
    return Lum/L_sol


"""
Star class that holds various properties
"""
class star:
    def __init__(self, name, temp, rstar, plx):

        self.name = name

        #Temperature in Kelvin
        self.Temp = temp

        #Radius in Rsol
        self.Rstar = rstar

        #parallax in mas
        self.plx = plx

        #dist in parsec
        self.dist = 1000/self.plx

        #Flux density in given wavelength filter (set in parameters above)
        self.flux = flux_density(temp,rstar,self.dist,wavelength,d_wavelength)

        #Bolometric Luminosity in Lsol
        self.Lstar = Luminosity(temp,rstar)

        #Habitable zone angle in mas
        self.hz_angle = (self.Lstar**.5)*plx


def get_nuller_response(M, telescope_positions, extent=3):

    angles = np.linspace(-extent*np.pi,extent*np.pi,sz)
    xy = np.meshgrid(angles, angles, indexing='ij')

    x = telescope_positions[:,0]
    y = telescope_positions[:,1]

    #Response is the 5 output electric fields as a function of the position on the sky
    response = np.zeros((M.shape[0],sz,sz), dtype='complex')

    for i in range(M.shape[1]):
        for k in range(M.shape[0]):
            #Inputs have a phase equal to xy array - linear multiplied by spatial frequency
            response[k] += np.exp(1j*(xy[0]*x[i] + xy[1]*y[i]))*M[k,i] #

    return np.abs(response)**2 #return intensity

def symmetrical_polygon(N=5,scale=1):
    angles = np.linspace(0,2*np.pi,N+1)
    xs = baseline*np.sin(angles)
    ys = baseline*np.cos(angles)
    return np.array([xs,ys]).T[:-1]

def plot_map(k,figure_number,extent):
    plt.figure(figure_number)
    plt.clf()
    plt.imshow(k, extent=extent)
    plt.xlabel('Offset (mas)')
    plt.ylabel('Offset (mas)')
    plt.colorbar()
    plt.tight_layout()

def limb_darkening(r):
    mu = np.sqrt(1-r**2)
    return 1-0.47*(1-mu)-0.23*(1-mu)**2


############################################################################
#Equation 2 of 2018 paper
N = knull.make_nuller_mat5()
tau_cet = star("Tau Ceti", temp=5344, rstar=0.793, plx=273.96)

bl = 50 #Baseline - to make habitable zone at maximum in nuller (phase = pi)

response_size = 3 #Need better name - number of phase cycles to

#Used in plotting - converts phase to angle on sky
pix2mas = wavelength/bl*rad2mas*response_size/sz

telescope_array = symmetrical_polygon(N=5,scale=1)

response = get_nuller_response(N, telescope_array, extent=response_size)

response /= (np.max(response[0])/5) #Response 0 is star output. So this is normalised by flux per telescope

#Create the kernel nulls. This is turning the output intensities into the kernel nulls (K in 2018 paper)
k1 = response[2]-response[3]
k2 = response[1]-response[4]


#Find the null function numerically
#Integral over theta, divide by two pi, as a function of radius
r_ix,y2 = azimuthalAverage(response[1], returnradii=True, binsize=0.8)
r_ix,y4 = azimuthalAverage(response[2], returnradii=True, binsize=0.8)
#averaging response function over theta as a function of radius. yi is average value


#Dodgy fit for a quartic and parabola
second_order_coeff = np.median(y2[1:16]/r_ix[1:16]**2)/mas_pix**2
fourth_order_coeff = np.median(y4[1:16]/r_ix[1:16]**2)/mas_pix**2


r = np.linspace(0,1,200)

#Total leakage is an integral of coeff multiplied by r^2 or r^4 and limb darkening intensity
mn_r2 = (np.trapz(I*r**3, r)/np.trapz(I*r, r))**.5
mn_r4 = (np.trapz(I*r**5, r)/np.trapz(I*r, r))**.25

#tau Ceti calculation (Scaling by radius of star)
star_mas = Rstar/214*plx
leakage_2nd = second_order_coeff*(mn_r2*star_mas)**2
leakage_4th = fourth_order_coeff*(mn_r4*star_mas)**4

#Convert from r to one quadrant of an xy grid
xy = np.meshgrid(r + 0.5*(r[1]-r[0]),r+ 0.5*(r[1]-r[0]))
rr = np.sqrt(xy[0]**2+xy[1]**2)
leakage_bracewell_const = 0.5*np.pi**2*np.sum(np.interp(rr,r,I)*xy[0]**2)/np.sum(np.interp(rr,r,I)) #Thing at top over x and y
#does interp go to 0 beyond r?
leakage_bracewell = leakage_bracewell_const*(star_mas/hz_angle)**2

#background in mag/sr
zodiacal_background =sfdsdfSD

#Leakage/zodiacal are as a fraction of the starlight that gets into the nulled output



#############################################################################

#Planet signal flux (photons/m^2/s)
signal


#Stellar leakage flux (photons/m^2/s)
leakage


#Zodiacal flux (photons/m^2/s)
zodiacal



################################################################################
C_per_t = telescope_area*exposure_time
#Signal
signal_per_t = np.sum(signal*C_per_t)
#Noise contributions: shot, leakage and zodiacal background
noise_per_t = np.sqrt(np.sum(signal*C_per_t) + np.sum(leakage*C_per_t) + zodiacal*C_per_t)
SNR_per_t = signal/noise


C_tot = total_area*exposure_time
#Signal
signal_tot = np.sum(signal*C_tot)
#Noise contributions: shot, leakage and zodiacal background
noise_tot = np.sqrt(np.sum(signal*C_tot) + np.sum(leakage*C_tot) + zodiacal*C_tot)
SNR_tot = signal/noise

###############################################################################

print("{:.1e} {:.1e} {:.1e} {:.1e}".format(leakage_2nd[0], leakage_4th[0], leakage_bracewell[0], zodiacal[0]))

#Stellar leakage domination over zodiacal
plt.figure(4)
plt.clf()
plt.semilogy(dist, leakage_2nd/zodiacal, color + '-', label=label + ' 5T Second Order')
plt.semilogy(dist, leakage_bracewell/zodiacal, color + '--', label=label + ' Bracewell')
plt.xlabel('Distance (pc)')
plt.ylabel('Leakage (fraction of Zodiacal)')


plt.figure(5)
#Variance in background divided by square of planet flux
var1 = (leakage_2nd + zodiacal)/(0.5*0.4**2) #Total noise flux, two fifths of planet light on average over star positions goes into nulled outputs, 0.5 comes from doubling the noise from the difference
var2 = (leakage_4th + zodiacal)/(0.5*0.4**2)
y = np.sqrt(zodiacal[0]*(1/var1 + 1/var2))
ref_dist = 1000./ref_plx

plt.plot(ref_dist, np.interp(ref_dist, dist[::-1], y[::-1]),color+'o')
plt.text(ref_dist+1, np.interp(ref_dist, dist[::-1], y[::-1]), star)
varz = (leakage_bracewell + zodiacal)/0.35**2
plt.plot(dist, y, color + '-', label=label + ' 5T Combined')
plt.plot(dist, np.sqrt(zodiacal[0]/varz), color + '--', label=label + ' Bracewell')
plt.ylabel('Relative SNR (equal total area)') #signal to noise compared to full zodiacal limit
plt.xlabel('Distance (pc)')
if not initial_run:
    plt.legend()
    plt.tight_layout()
