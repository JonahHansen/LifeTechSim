import sys
sys.path.append("../..")
import numpy as np
from opticstools import knull

#The baseline is the distance between adjacent spacecraft
def triangle(baseline):
    R = 0.5773*baseline
    angles = np.linspace(0,2*np.pi,4)
    xs = R*np.cos(angles)
    ys = R*np.sin(angles)
    return np.array([xs,ys]).T[:-1]


def get_nuller_response(baseline,fov,sz,base_wavelength):

    N = knull.make_nuller_mat3()

    telescope_array = triangle(baseline)

    sky_angles = np.linspace(-fov/2,fov/2,sz)

    xy = np.meshgrid(sky_angles, sky_angles, indexing='ij')

    #x,y are telescope positions in units of wavelength
    x = telescope_array[:,0]/base_wavelength
    y = telescope_array[:,1]/base_wavelength

    #Response is the 5 output electric fields as a function of the position on the sky
    response = np.zeros((3,sz,sz), dtype='complex')

    for i in range(3):
        for k in range(3):
            #Inputs have a phase equal to xy array - linear multiplied by spatial frequency
            #ul + vm?
            response[k] += np.exp(2*np.pi*1j*(xy[0]*x[i] + xy[1]*y[i]))*N[k,i] #

    response = np.abs(response)**2
    response /= (np.max(response[0])/3) #normalise by flux per telescope

    #Create the kernel nulls. This is turning the output intensities into the kernel nulls (K in 2018 paper)
    k = response[1]-response[2]

    return [(response[1], k)] #return intensity per telescope
