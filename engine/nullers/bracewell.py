import sys
sys.path.append("../..")
import numpy as np

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
Function to calculate the response map of a four telescope bracewell-type interferometer.

Inputs:
    baseline = length of the shorter, nulling baseline in meters
    fov = total field of view of the interferometer in radians
    sz = size of the array
    base_wavelength = wavelength upon which the baseline is optimised. In meters
    ratio = ratio of the imaging to nulling baseline. Default is 6

Outputs:
    List of modulation maps of the form:
        (Transmission map, Kernel map)
    The transmission map is one of the nulled outputs that is chosen to create the kernel.
"""
def get_nuller_response(baseline,fov,sz,base_wavelength,ratio=6):

    #Two couplers that combine the nulling baseline telescopes
    R1 = 1/np.sqrt(2)*np.array([[1,1,0,0],
                                [1,-1,0,0],
                                [0,0,1,-1],
                                [0,0,1,1]])

    #Phase chopping the two pairs of nulls
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
            #ul + vm
            response[k] += np.exp(2*np.pi*1j*(xy[0]*x[i] + xy[1]*y[i]))*M[k,i] #

    response = np.abs(response)**2 #To intensity

    k = response[1] - response[2] #Kernel is the difference

    return [(response[1],k)] #return intensity per telescope
