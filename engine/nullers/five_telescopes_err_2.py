import numpy as np

"""
K-5 nuller, except with errors introduced (see paper 7)

"""

"""
Calculates the positions of five telescopes in a regular pentagonal formation.

Inputs:
    baseline = the shorter baseline (distance between adjacent telescopes)
Output:
    List of x positions and list of y positions
"""
def pentagon(baseline):
    R = 0.85065*baseline #To the radius of the pentagon
    angles = np.linspace(0,2*np.pi,6) #Angular positions of the telescopes
    #To cartesian
    xs = R*np.cos(angles)
    ys = R*np.sin(angles)
    return np.array([xs,ys]).T[:-1]


"""
Calculate the beam splitting module matrix

Inputs:
    n = size of beam combiner matrix
    i = location of the beam splitter in the combiner
    phi = phase shifting angle
    theta = mixing angle for beamsplitter
"""
def calc_a(n,i,phi,theta):
    a = np.identity(n,complex)

    a[i,i] = np.sin(theta)
    a[i,i+1] = np.exp(1j*phi)*np.cos(theta)
    a[i+1,i] = np.cos(theta)
    a[i+1,i+1] = -np.exp(1j*phi)*np.sin(theta)

    return a



"""
Calculate the beam combiner matrix

Inputs:
    dphi = error in phase shifting angle
    dR = error in beam splitter reflectance
"""
def calc_K5_M(dphi,dR):

    #Mixing angles
    theta = np.array([7*np.pi/4,
                     np.arcsin(1/np.sqrt(3)),
                     np.arctan(np.sqrt(1/3*(4-np.sqrt(5)))),
                     5*np.pi/6,
                     np.pi-np.arcsin(1/3*np.sqrt(3-2/np.sqrt(5))),
                     np.pi+np.arcsin((-1+3*np.sqrt(5))/(2*np.sqrt(22))),
                     np.pi-np.arcsin(1/np.sqrt(5)),
                     -np.pi/6,
                     -np.arcsin(1/np.sqrt(3)),
                     -np.pi/4])

    #Phase shifting angles
    phi = np.array([np.pi,
                     np.pi,
                     -3*np.pi/10-np.arctan(np.sqrt(5-2/np.sqrt(5))),
                     np.pi,
                     -np.pi+np.arctan(3/7*np.sqrt(5-2*np.sqrt(5))),
                     np.pi-np.arctan(1/(np.sqrt(5*(5+2*np.sqrt(5))))),
                     np.pi,
                     np.arctan(np.sqrt(2-2/np.sqrt(5))),
                     -np.arctan(np.sqrt(1/10*(5-np.sqrt(5)))),
                     -np.pi/2])

    #Convert reflectance error to mixing angle error
    dtheta = dR/np.cos(theta)

    #Beam splitter matrices
    a1 = calc_a(5,3,phi[0]+dphi[0],theta[0]+dtheta[0])
    a2 = calc_a(5,2,phi[1]+dphi[1],theta[1]+dtheta[1])
    a3 = calc_a(5,3,phi[2]+dphi[2],theta[2]+dtheta[2])
    a4 = calc_a(5,1,phi[3]+dphi[3],theta[3]+dtheta[3])
    a5 = calc_a(5,2,phi[4]+dphi[4],theta[4]+dtheta[4])
    a6 = calc_a(5,3,phi[5]+dphi[5],theta[5]+dtheta[5])
    a7 = calc_a(5,0,phi[6]+dphi[6],theta[6]+dtheta[6])
    a8 = calc_a(5,1,phi[7]+dphi[7],theta[7]+dtheta[7])
    a9 = calc_a(5,2,phi[8]+dphi[8],theta[8]+dtheta[8])
    a10 = calc_a(5,3,phi[9]+dphi[9],theta[9]+dtheta[9])

    #Phase shift at output
    C = np.identity(5,complex)
    C[4,4] = np.exp(1j*np.pi)

    ls_mats = [C,a10,a9,a8,a7,a6,a5,a4,a3,a2,a1]

    M = np.identity(5)
    for ai in ls_mats:
        M = np.matmul(M,ai)

    return M

"""
Function to calculate the response map of a five telescope kernel nuller interferometer with errors (see paper 7)

Inputs:
    dphi = error in phase shift plate of the beam combiner
    dR = error in the reflectance of a beam splitter
    baseline = length of the shorter baseline in meters (adjacent telescopes)
    fov = total field of view of the interferometer in radians
    sz = size of the array
    base_wavelength = wavelength upon which the baseline is optimised. In meters

Outputs:
    List of modulation maps of the form:
        (Transmission map, Kernel map)
    The transmission map is one of the nulled outputs that is chosen to create the kernel.
"""
def get_nuller_response(dphi,dR,baseline,fov,sz,base_wavelength):

    M = calc_K5_M(dphi,dR)

    telescope_array = pentagon(baseline)

    sky_angles = np.linspace(-fov/2,fov/2,sz)

    xy = np.meshgrid(sky_angles, sky_angles, indexing='ij')

    #x,y are telescope positions in units of wavelength
    x = telescope_array[:,0]/base_wavelength
    y = telescope_array[:,1]/base_wavelength

    #Response is the 5 output electric fields as a function of the position on the sky
    response = np.zeros((5,sz,sz), dtype='complex')

    for i in range(5):
        for k in range(5):
            #Inputs have a phase equal to xy array - linear multiplied by spatial frequency
            #ul + vm
            response[k] += np.exp(2*np.pi*1j*(xy[0]*x[i] + xy[1]*y[i]))*M[k,i] #

    response = np.abs(response)**2 #To intensity
    response /= (np.max(response[0])/5) #normalise by flux per telescope

    #This is turning the output intensities into the kernel nulls (K in 2018 paper)
    k1 = response[3]-response[4]
    k2 = response[1]-response[2]

    return [(response[1],k1),(response[2],k2)] #return intensity per telescope
