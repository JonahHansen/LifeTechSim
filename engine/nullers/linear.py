import numpy as np



def linear(baseline_short):
    xs = np.array([-2,-1,1,2])*baseline_short
    ys = [0,0,0,0]
    return np.array([xs,ys]).T


def guyon_fig3_mat4_redone():

    M = np.array([[0.5, 0.5, 0.5, 0.5],
                 [np.sqrt(2/5), 1/(np.sqrt(10)), -1/(np.sqrt(10)), -np.sqrt(2/5)],
                 [np.sqrt(7/40)*np.exp(1j*np.pi), np.sqrt(13/40)*np.exp(0.107554*1j*np.pi), np.sqrt(13/40)*np.exp(-0.466571*1j*np.pi), np.sqrt(7/40)*np.exp(0.640983*1j*np.pi)],
                 [np.sqrt(7/40)*np.exp(1j*np.pi), np.sqrt(13/40)*np.exp(-0.107554*1j*np.pi), np.sqrt(13/40)*np.exp(0.466571*1j*np.pi), np.sqrt(7/40)*np.exp(-0.640983*1j*np.pi)]], dtype=complex)

    return M


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
    M = guyon_fig3_mat4_redone()

    telescope_array = linear(baseline)

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
            #ul + vm
            response[k] += np.exp(2*np.pi*1j*(xy[0]*x[i] + xy[1]*y[i]))*M[k,i] #

    response = np.abs(response)**2 #To intensities
    response /= (np.max(response[0])/4) #normalise by flux per telescope

    #Create the kernel nulls. This is turning the output intensities into the kernel nulls (K in 2018 paper)
    k1 = response[2]-response[3]

    return [(response[1],k1)] #return intensity per telescope
