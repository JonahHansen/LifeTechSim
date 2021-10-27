import sys
sys.path.append("../..")
import numpy as np

#The baseline given to this function is the NULLING baseline
#The ratio is how many times longer the imaging baseline is
def rectangle(baseline,ratio):
    xs = [baseline/2,-baseline/2,-baseline/2,baseline/2]
    ys = [baseline*ratio/2,baseline*ratio/2,-baseline*ratio/2,-baseline*ratio/2]
    return np.array([xs,ys]).T


def get_nuller_response(baseline,fov,sz,base_wavelength,ratio=6):

    R1 = 1/np.sqrt(2)*np.array([[1,1,0,0],
                                [1,-1,0,0],
                                [0,0,1,-1],
                                [0,0,1,1]])

    R2 = 1/np.sqrt(2)*np.array([[np.sqrt(2),0,0,0],
                                [0,1,1j,0],
                                [0,1j,1,0],
                                [0,0,0,np.sqrt(2)]])

    M = np.matmul(R2,R1)

    telescope_array = rectangle(baseline,ratio)

    sky_angles = np.linspace(-fov/2,fov/2,sz)

    xy = np.meshgrid(sky_angles, sky_angles, indexing='ij')

    #x,y are telescope positions in units of wavelength
    x = telescope_array[:,0]/base_wavelength
    y = telescope_array[:,1]/base_wavelength

    #Response is the 5 output electric fields as a function of the position on the sky
    response = np.zeros((4,sz,sz), dtype='complex')

    for i in range(4):
        for k in range(4):
            #Inputs have a phase equal to xy array - linear multiplied by spatial frequency
            #ul + vm?
            response[k] += np.exp(2*np.pi*1j*(xy[0]*x[i] + xy[1]*y[i]))*M[k,i] #

    response = np.abs(response)**2

    k = response[1] - response[2]

    return [(response[1],k)] #return intensity per telescope
