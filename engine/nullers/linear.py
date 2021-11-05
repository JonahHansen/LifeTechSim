import sys
sys.path.append("../..")
import numpy as np
from opticstools import knull

"""
Calculates the positions of four telescopes in a rectangular formation.

Inputs:
    baseline = the shorter baseline (ie nulling baseline)
    ratio = how many times longer the imaging baseline is
Output:
    List of x positions and list of y positions
"""
def linear(baseline_short):
    xs = np.array([-3,3,1,-2])*baseline_short
    ys = [0,0,0,0]
    return np.array([xs,ys]).T


def guyon_fig6_mat4(phase_shifters=np.array([1,1,1], dtype=complex), ker_coupler_angle=np.pi/2):
    """Make a nuller electric field matrix for a 4 telescope combiner

    Based on Guyon 2013, with kernel nulling mixing matrix

    phase_shifters: array-like
        Three phasors to represent the on-chip phase modulation. For no
        phase shift, these should be 1, and for 180 degrees they should
        be -1 (i.e. interferometric chopping)

    ker_coupler_angle: float
        In the case of kernel outputs only, the final 50/50 couplers are defined by the
        one phase. See e.g. https://en.wikipedia.org/wiki/Unitary_matrix, setting theta
        to pi/2 (50/50), and phase angle varphi_1 to 0.
    """

    MM0 = np.array([[np.sqrt(0.25), np.sqrt(0.25), np.sqrt(0.25), np.sqrt(0.25)],
                    [np.sqrt(0.3324), -np.sqrt(0.4643), np.sqrt(0.0687), -np.sqrt(0.1346)],
                    [np.sqrt(0.2057), -np.sqrt(0.2518), -np.sqrt(0.4695), np.sqrt(0.0730)],
                    [np.sqrt(0.2119), np.sqrt(0.0339), -np.sqrt(0.2119), -np.sqrt(0.5424)]], dtype=complex)

    #Add in the phase shifters
    for ix,phasor in enumerate(phase_shifters):
        MM0[ix+1] *= phasor

    MM1 = np.zeros( (7,4), dtype=complex)
    #Start with the bright output.
    MM1[0] = MM0[0]

    #Now lets take the three nulled outputs, and put each of these into a 2x2-coupler.
    PHI0 = np.exp(1j*ker_coupler_angle)
    PHI1 = np.conj(PHI0)

    MM1[1] = (MM0[1] + PHI0*MM0[2]) / 2
    MM1[2] = (-PHI1*MM0[1] + MM0[2]) / 2

    MM1[3] = (MM0[1] + PHI0*MM0[3]) / 2
    MM1[4] = (-PHI1*MM0[1] + MM0[3]) / 2

    MM1[5] = (MM0[2] + PHI0*MM0[3]) / 2
    MM1[6] = (-PHI1*MM0[2] + MM0[3]) / 2

    return MM1


"""
Function to calculate the response map of a four telescope kernel nuller interferometer.

Inputs:
    baseline = length of the shorter baseline in meters
    fov = total field of view of the interferometer in radians
    sz = size of the array
    base_wavelength = wavelength upon which the baseline is optimised. In meters

Outputs:
    List of modulation maps of the form:
        (Transmission map, Kernel map)
    The transmission map is one of the nulled outputs that is chosen to create the kernel.
"""
def get_nuller_response(baseline,fov,sz,base_wavelength):

    #Nulling matrix from Martinache and Ireland (2008) paper
    M = guyon_fig6_mat4()

    telescope_array = linear(baseline)

    sky_angles = np.linspace(-fov/2,fov/2,sz)

    xy = np.meshgrid(sky_angles, sky_angles, indexing='ij')

    #x,y are telescope positions in units of wavelength
    x = telescope_array[:,0]/base_wavelength
    y = telescope_array[:,1]/base_wavelength

    #Response is the 5 output electric fields as a function of the position on the sky
    response = np.zeros((7,sz,sz), dtype='complex')

    for i in range(4):
        for k in range(7):
            #Inputs have a phase equal to xy array - linear multiplied by spatial frequency
            #ul + vm
            response[k] += np.exp(2*np.pi*1j*(xy[0]*x[i] + xy[1]*y[i]))*M[k,i] #

    response = np.abs(response)**2 #To intensities
    response /= (np.max(response[0])/4) #normalise by flux per telescope

    #Create the kernel nulls. This is turning the output intensities into the kernel nulls (K in 2018 paper)
    k1 = response[1]-response[2]
    k2 = response[3]-response[4]
    k3 = response[5]-response[6]

    return [(response[1],k1),(response[3],k2),(response[5],k3)] #return intensity per telescope
