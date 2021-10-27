import sys
sys.path.append("../..")
import numpy as np
from opticstools import knull

def rectangle(baseline,ratio):
    xs = [baseline/2,-baseline/2,-baseline/2,baseline/2]
    ys = [baseline*ratio/2,baseline*ratio/2,-baseline*ratio/2,-baseline*ratio/2]
    return np.array([xs,ys]).T

#The baseline given to this function is the one that defines the side length of the pentagon
#Scale the argument to this function appropriately if defining baselines based on diagonals
def right_kite(baseline,ratio):

    R = baseline*np.sqrt(1+ratio**2)/2
    circ_angle = 2*np.arctan(1/ratio)
    angles = np.array([np.pi/2-circ_angle,np.pi/2,np.pi/2+circ_angle,3*np.pi/2])
    xs = R*np.cos(angles)
    ys = R*np.sin(angles)
    return np.array([xs,ys]).T


def get_nuller_response(baseline,fov,sz,base_wavelength,ratio=1.69):

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
            #ul + vm?
            response[k] += np.exp(2*np.pi*1j*(xy[0]*x[i] + xy[1]*y[i]))*M[k,i] #

    response = np.abs(response)**2
    response /= (np.max(response[0])/4) #normalise by flux per telescope

    #Create the kernel nulls. This is turning the output intensities into the kernel nulls (K in 2018 paper)

    k1 = response[1]-response[2]
    k2 = response[3]-response[4]
    k3 = response[5]-response[6]

    return [(response[1],k1),(response[3],k2),(response[5],k3)] #return intensity per telescope
