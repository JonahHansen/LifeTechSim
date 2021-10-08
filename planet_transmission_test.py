import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as const
from opticstools import knull

baseline = 80 #m
wavelength = 10e-6 #m
scale_factor = 3
sz = 400

rad2mas = np.degrees(1)*3600e3 #Number of milliarcsec in one radian

#Tau Ceti
L = 0.52 #Lsol
Dist = 3.65 #Pc

HZAngle = np.sqrt(L)*1000/Dist
baseline = wavelength/2*rad2mas/HZAngle
pix2mas = wavelength/baseline*rad2mas*scale_factor/sz


def azimuthal_rms(image,r):
    n_angles = 10000
    angles = np.linspace(0,2*np.pi,n_angles)

    centre = (int(image.shape[0]/2),int(image.shape[1]/2))
    sum = 0
    for theta in angles:
        x = centre[0] + r*np.cos(theta)
        y = centre[1] + r*np.sin(theta)

        a = image[int(x),int(y)]

        sum += a**2

    return np.sqrt(sum/n_angles)


def pentagon(baseline):
    angles = np.linspace(0,2*np.pi,6)
    xs = 1.176*baseline*np.sin(angles)
    ys = 1.176*baseline*np.cos(angles)
    return np.array([xs,ys]).T[:-1]

def get_nuller_response():

    N = knull.make_nuller_mat5()

    telescope_array = pentagon(baseline)

    fov = wavelength/baseline*scale_factor #?????
    sky_angles = np.linspace(-fov/2,fov/2,sz)

    xy = np.meshgrid(sky_angles, sky_angles, indexing='ij')

    #x,y are telescope positions in units of wavelength
    x = telescope_array[:,0]/wavelength
    y = telescope_array[:,1]/wavelength

    #Response is the 5 output electric fields as a function of the position on the sky
    response = np.zeros((5,sz,sz), dtype='complex')

    for i in range(5):
        for k in range(5):
            #Inputs have a phase equal to xy array - linear multiplied by spatial frequency
            #ul + vm?
            response[k] += np.exp(2*np.pi*1j*(xy[0]*x[i] + xy[1]*y[i]))*N[k,i] #

    response = np.abs(response)**2
    response /= (np.max(response[0])/5) #Normalise by flux per telescope

    #Create the kernel nulls. This is turning the output intensities into the kernel nulls (K in 2018 paper)

    k1 = response[2]-response[3]
    k2 = response[1]-response[4]

    return response, k1, k2 #return intensity per telescope

r,k1,k2 = get_nuller_response()

rs = np.linspace(0.01,199)
k1_ave = []
k2_ave = []

for r in rs:
    k1_ave.append(azimuthal_rms(k1,r))
    k2_ave.append(azimuthal_rms(k2,r))

plt.figure(1)
plt.imshow(k1)
plt.figure(2)
plt.imshow(k2)
plt.figure(3)
plt.plot(rs*pix2mas,k1_ave,label="k1")
plt.plot(rs*pix2mas,k2_ave,label="k2")
plt.axvline(x=HZAngle,color="k")
plt.legend()
plt.xlabel("Radial coordinate (mas)")
plt.ylabel("Transmission")
