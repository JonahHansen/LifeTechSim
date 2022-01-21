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

def triangle(baseline):
    R = 0.5773*baseline #Convert to circumcircle radius
    angles = np.linspace(0,2*np.pi,4) #angular coordinates of the telescopes
    xs = R*np.cos(angles)
    ys = R*np.sin(angles)
    return np.array([xs,ys]).T[:-1]

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
def get_nuller_response_53(baseline,fov,sz,base_wavelength):

    a = 0.9048040910575242
    b = 0.7551280988643292
    c = 0.2707664135274219
    d = 1.1469849337259779
    e = 0.6340376775301025
    f = 0.3918568348616487

    N = 1/np.sqrt(5)*np.array([[1,1,1],
                                [-np.sqrt(5)/6+a*1j,-np.sqrt(5)/6-c*1j,np.sqrt(5)/3-e*1j],
                                [-np.sqrt(5)/6+b*1j,-np.sqrt(5)/6-d*1j,np.sqrt(5)/3+f*1j],
                                [-np.sqrt(5)/6-b*1j,-np.sqrt(5)/6+d*1j,np.sqrt(5)/3-f*1j],
                                [-np.sqrt(5)/6-a*1j,-np.sqrt(5)/6+c*1j,np.sqrt(5)/3+e*1j]])


    telescope_array = triangle(baseline)

    sky_angles = np.linspace(-fov/2,fov/2,sz)

    xy = np.meshgrid(sky_angles, sky_angles, indexing='ij')

    #x,y are telescope positions in units of wavelength
    x = telescope_array[:,0]/base_wavelength
    y = telescope_array[:,1]/base_wavelength

    #Response is the 5 output electric fields as a function of the position on the sky
    response = np.zeros((5,sz,sz), dtype='complex')

    for i in range(3):
        for k in range(5):
            #Inputs have a phase equal to xy array - linear multiplied by spatial frequency
            #ul + vm
            response[k] += np.exp(2*np.pi*1j*(xy[0]*x[i] + xy[1]*y[i]))*N[k,i] #

    response = np.abs(response)**2 #To intensity
    response /= (np.max(response[0])/5) #normalise by flux per telescope

    #This is turning the output intensities into the kernel nulls (K in 2018 paper)
    k1 = response[1]-response[4]
    k2 = response[2]-response[3]

    return [(response[1],k1),(response[2],k2)] #return intensity per telescope




def get_nuller_response_54(baseline,fov,sz,base_wavelength,ratio=6):

    a = 0.9510565162951535
    b = 0.5877852522924732

    N = 1/np.sqrt(5)*np.array([[1,1,1,1],
                                [np.sqrt(5)/4+a*1j,-np.sqrt(5)/4+b*1j,-np.sqrt(5)/4-b*1j,np.sqrt(5)/4-a*1j],
                                [np.sqrt(5)/4-b*1j,-np.sqrt(5)/4+a*1j,-np.sqrt(5)/4-a*1j,np.sqrt(5)/4+b*1j],
                                [np.sqrt(5)/4+b*1j,-np.sqrt(5)/4-a*1j,-np.sqrt(5)/4+a*1j,np.sqrt(5)/4-b*1j],
                                [np.sqrt(5)/4-a*1j,-np.sqrt(5)/4-b*1j,-np.sqrt(5)/4+b*1j,np.sqrt(5)/4+a*1j]])

    telescope_array = rectangle(baseline,ratio)

    sky_angles = np.linspace(-fov/2,fov/2,sz)

    xy = np.meshgrid(sky_angles, sky_angles, indexing='ij')

    #x,y are telescope positions in units of wavelength
    x = telescope_array[:,0]/base_wavelength
    y = telescope_array[:,1]/base_wavelength

    #Response is the 5 output electric fields as a function of the position on the sky
    response = np.zeros((5,sz,sz), dtype='complex')

    for i in range(4):
        for k in range(5):
            #Inputs have a phase equal to xy array - linear multiplied by spatial frequency
            #ul + vm
            response[k] += np.exp(2*np.pi*1j*(xy[0]*x[i] + xy[1]*y[i]))*N[k,i] #

    response = np.abs(response)**2 #To intensity
    response /= (np.max(response[0])/5) #normalise by flux per telescope

    #This is turning the output intensities into the kernel nulls (K in 2018 paper)
    k1 = response[1]-response[4]
    k2 = response[2]-response[3]

    return [(response[1],k1),(response[2],k2)] #return intensity per telescope


def get_nuller_response_43(baseline,fov,sz,base_wavelength):

    M = 1/np.sqrt(4)*np.array([[1,1,1],
                                [1/3+1j,1/3-1j,-2/3],
                                [1/3-1j,1/3+1j,-2/3],
                                [2/3,2/3,-4/3]])

    telescope_array = triangle(baseline)

    sky_angles = np.linspace(-fov/2,fov/2,sz)

    xy = np.meshgrid(sky_angles, sky_angles, indexing='ij')

    #x,y are telescope positions in units of wavelength
    x = telescope_array[:,0]/base_wavelength
    y = telescope_array[:,1]/base_wavelength

    #Response is the 5 output electric fields as a function of the position on the sky
    response = np.zeros((4,sz,sz), dtype='complex')

    for i in range(3):
        for k in range(4):
            #Inputs have a phase equal to xy array - linear multiplied by spatial frequency
            #ul + vm
            response[k] += np.exp(2*np.pi*1j*(xy[0]*x[i] + xy[1]*y[i]))*M[k,i] #

    response = np.abs(response)**2 #To intensity

    k = response[1] - response[2] #Kernel is the difference

    return [(response[1],k)] #return intensity per telescope
