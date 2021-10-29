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
def rectangle(baseline,ratio):
    xs = [baseline/2,-baseline/2,-baseline/2,baseline/2]
    ys = [baseline*ratio/2,baseline*ratio/2,-baseline*ratio/2,-baseline*ratio/2]
    return np.array([xs,ys]).T

"""
Calculates the positions of four telescopes in a right kite formation.
See PDF notes for why this works

Inputs:
    baseline = the shorter baseline (ie top half of the kite)
    ratio = how many times longer the bottom sides of the kite are
Output:
    List of x positions and list of y positions
"""
def right_kite(baseline,ratio):

    R = baseline*np.sqrt(1+ratio**2)/2 #radius of the telescope positions
    circ_angle = 2*np.arctan(1/ratio) #Angle between the top telescope and either of the side telescopes
    angles = np.array([np.pi/2-circ_angle,np.pi/2,np.pi/2+circ_angle,3*np.pi/2]) #Angular positions
    xs = R*np.cos(angles)
    ys = R*np.sin(angles)
    return np.array([xs,ys]).T

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
def get_nuller_response(baseline,fov,sz,base_wavelength,ratio=1.69):

    #Nulling matrix from Martinache and Ireland (2008) paper
    M = knull.make_nuller_mat4(bracewell_design=False, ker_only=True)

    telescope_array = right_kite(baseline,ratio)

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
