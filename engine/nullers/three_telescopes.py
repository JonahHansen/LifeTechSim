import numpy as np

"""
Calculates the positions of five telescopes in an equilateral triangle formation.

Inputs:
    baseline = the distance between adjacent telescopes
Output:
    List of x positions and list of y positions
"""
def triangle(baseline):
    R = 0.5773*baseline #Convert to circumcircle radius
    angles = np.linspace(0,2*np.pi,4) #angular coordinates of the telescopes
    xs = R*np.cos(angles)
    ys = R*np.sin(angles)
    return np.array([xs,ys]).T[:-1]

def make_nuller_mat3_JH():
    """Create a 3x3 Nuller matrix"""
    initial_mat = np.zeros( (3,3) )
    for i in range(3):
        initial_mat[i] = np.arange(3) * i
    initial_mat = np.exp(2j*np.pi/3*initial_mat)
    return initial_mat


"""
Function to calculate the response map of a three telescope kernel nuller interferometer.

Inputs:
    baseline = length of the baseline in meters
    fov = total field of view of the interferometer in radians
    sz = size of the array
    base_wavelength = wavelength upon which the baseline is optimised. In meters

Outputs:
    List of modulation maps of the form:
        (Transmission map, Kernel map)
    The transmission map is one of the nulled outputs that is chosen to create the kernel.
"""
def get_nuller_response(baseline,fov,sz,base_wavelength):

    N = make_nuller_mat3_JH()

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
            #ul + vm
            response[k] += np.exp(2*np.pi*1j*(xy[0]*x[i] + xy[1]*y[i]))*N[k,i] #

    response = np.abs(response)**2 #To intensity

    response /= (np.max(response[0])/3) #normalise by flux per telescope

    #Create the kernel nulls. This is turning the output intensities into the kernel nulls (K in 2018 paper)
    k = response[1]-response[2]

    return [(response[1], k)] #return intensity per telescope
