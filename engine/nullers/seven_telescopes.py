import numpy as np

"""
Calculates the positions of seven telescopes in a regular heptagonal formation.

Inputs:
    baseline = the shorter baseline (distance between adjacent telescopes)
Output:
    List of x positions and list of y positions
"""
def heptagon(baseline):
    R = 1.152*baseline #To the radius of the pentagon
    angles = np.linspace(0,2*np.pi,8) #Angular positions of the telescopes
    #To cartesian
    xs = R*np.cos(angles)
    ys = R*np.sin(angles)
    return np.array([xs,ys]).T[:-1]

def make_nuller_mat7():
    """Create a 7x7 Nuller matrix"""
    initial_mat = np.zeros( (7,7) )
    for i in range(7):
        initial_mat[i] = np.arange(7) * i
    initial_mat = np.exp(2j*np.pi/7*initial_mat)
    return initial_mat

"""
Function to calculate the response map of a seven telescope kernel nuller interferometer.

Inputs:
    baseline = length of the shorter baseline in meters (adjacent telescopes)
    fov = total field of view of the interferometer in radians
    sz = size of the array
    base_wavelength = wavelength upon which the baseline is optimised. In meters

Outputs:
    List of modulation maps of the form:
        (Transmission map, Kernel map)
    The transmission map is one of the nulled outputs that is chosen to create the kernel.
"""
def get_nuller_response(baseline,fov,sz,base_wavelength):

    N = make_nuller_mat7()

    telescope_array = heptagon(baseline)

    sky_angles = np.linspace(-fov/2,fov/2,sz)

    xy = np.meshgrid(sky_angles, sky_angles, indexing='ij')

    #x,y are telescope positions in units of wavelength
    x = telescope_array[:,0]/base_wavelength
    y = telescope_array[:,1]/base_wavelength

    #Response is the 5 output electric fields as a function of the position on the sky
    response = np.zeros((7,sz,sz), dtype='complex')

    for i in range(7):
        for k in range(7):
            #Inputs have a phase equal to xy array - linear multiplied by spatial frequency
            #ul + vm
            response[k] += np.exp(2*np.pi*1j*(xy[0]*x[i] + xy[1]*y[i]))*N[k,i] #

    response = np.abs(response)**2 #To intensity
    response /= (np.max(response[0])/7) #normalise by flux per telescope

    #This is turning the output intensities into the kernel nulls (K in 2018 paper)
    k1 = response[1]-response[6]
    k2 = response[2]-response[5]
    k3 = response[3]-response[4]

    import matplotlib.pyplot as plt
    for i in range(7):
        plt.figure(i)
        plt.imshow(response[i])
    plt.show()


    return [(response[1],k1),(response[2],k2),(response[3],k3)] #return intensity per telescope
