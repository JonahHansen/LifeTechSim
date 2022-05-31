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

def make_nuller_mat4(bracewell_design=True, ker_only=False,
    phase_shifters=np.array([1,1,1], dtype=complex), ker_coupler_angle=np.pi/2):
    """Make a nuller electric field matrix for a 3 telescope combiner

    Parameters
    ----------
    bracewell_design: bool
        Do we use a Bracewell-like design? If not, then we use an
        obvious symmetrical matrix. Harry can build a bracewell-like
        design already.

    ker_only: bool
        Do we use kernel outputs only or all 9?

    phase_shifters: array-like
        Three phasors to represent the on-chip phase modulation. For no
        phase shift, these should be 1, and for 180 degrees they should
        be -1 (i.e. interferometric chopping)

    ker_coupler_angle: float
        In the case of kernel outputs only, the final 50/50 couplers are defined by the
        one phase. See e.g. https://en.wikipedia.org/wiki/Unitary_matrix, setting theta
        to pi/2 (50/50), and phase angle varphi_1 to 0.
    """
    #4x4 is nice and symmetric because you can get flux and visibilities
    #from 3 tri-couplers in the nulled channels.
    if bracewell_design:
        sq2 = np.sqrt(2)
        #For a "Bracewell" design, each element has a
        #matrix 1/sq2 * np.array([[1,1],[1,-1]])
        MM0 = 0.5 * np.array([[1,      1,   1,    1],
                              [1,      1,  -1,   -1],
                              [sq2, -sq2,   0,    0],
                              [0,      0, sq2, -sq2]], dtype=complex)
    else:
        MM0 = 0.5 * np.array([[1, 1, 1, 1],
                              [1, 1,-1,-1],
                              [1,-1, 1,-1],
                              [1,-1,-1, 1]], dtype=complex)

    #Add in the phase shifters
    for ix,phasor in enumerate(phase_shifters):
        MM0[ix+1] *= phasor

    if ker_only:
        MM1 = np.zeros( (7,4), dtype=complex)
        #Start with the bright output.
        MM1[0] = MM0[0]

        #Now lets take the three nulled outputs, and put each of these into a 2x2-coupler.
        PHI0 = np.exp(1j*ker_coupler_angle)
        PHI1 = np.conj(PHI0)

        MM1[1] = (MM0[1] + PHI0*MM0[2]) / 2
        MM1[2] = (-PHI1*MM0[1] + MM0[2]) / 2

        MM1[3] = (MM0[1] + PHI0*MM0[3]) / 2
        MM1[4] = (-PHI1*MM0[1] + MM0[3]) / 2

        MM1[5] = (MM0[2] + PHI0*MM0[3]) / 2
        MM1[6] = (-PHI1*MM0[2] + MM0[3]) / 2
    else:
        #Now lets take the three nulled outputs, and put each of these into a tri-coupler.
        PHI0 = np.exp(2j*np.pi/3)
        MM1 = np.zeros( (10,4), dtype=complex)
        MM1[0] = MM0[0]

        MM1[1] = (MM0[1] + MM0[2] * PHI0**0) / np.sqrt(6)
        MM1[2] = (MM0[1] + MM0[2] * PHI0**1) / np.sqrt(6)
        MM1[3] = (MM0[1] + MM0[2] * PHI0**2) / np.sqrt(6)

        MM1[4] = (MM0[1] + MM0[3] * PHI0**0) / np.sqrt(6)
        MM1[5] = (MM0[1] + MM0[3] * PHI0**1) / np.sqrt(6)
        MM1[6] = (MM0[1] + MM0[3] * PHI0**2) / np.sqrt(6)

        MM1[7] = (MM0[2] + MM0[3] * PHI0**0) / np.sqrt(6)
        MM1[8] = (MM0[2] + MM0[3] * PHI0**1) / np.sqrt(6)
        MM1[9] = (MM0[2] + MM0[3] * PHI0**2) / np.sqrt(6)
    return MM1


"""
Calculates the positions of four telescopes in a right kite formation.
See PDF notes for why this works

Inputs:
    baseline = the shorter baseline (ie top half of the kite)
    ratio = how many times longer the bottom sides of the kite are
Output:
    List of x positions and list of y positions
"""
def formation(baseline):
    R = 0.5773*baseline #Convert to circumcircle radius
    angles = np.linspace(0,2*np.pi,4) #angular coordinates of the telescopes
    xs = R*np.cos(angles)
    xs = np.insert(xs,0,0)
    ys = R*np.sin(angles)
    ys = np.insert(ys,0,0)
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
def get_nuller_response(baseline,fov,sz,base_wavelength):

    #Nulling matrix from Martinache and Ireland (2008) paper
    M = make_nuller_mat4(bracewell_design=False, ker_only=True)

    telescope_array = formation(baseline)

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
