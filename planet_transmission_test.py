import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as const
from opticstools import knull

baseline = 80 #m
wavelength = 14e-6 #m
scale_factor = 1.2
sz = 400

rad2mas = np.degrees(1)*3600e3 #Number of milliarcsec in one radian

#Tau Ceti
L = 0.2 #Lsol
Dist = 10.65 #Pc

HZAngle = np.sqrt(L)*1000/Dist
baseline = wavelength*rad2mas/HZAngle
fov = 2*scale_factor*HZAngle/rad2mas
pix2mas = fov*rad2mas/sz


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
    xs = 1/1.90211*2.08*baseline*np.sin(angles)
    ys = 1/1.90211*2.08*baseline*np.cos(angles)
    return np.array([xs,ys]).T[:-1]

def pentagon2(baseline):
    angles = np.linspace(0,2*np.pi,6)
    xs = 1/1.17557*baseline*np.sin(angles)
    ys = 1/1.17557*baseline*np.cos(angles)
    return np.array([xs,ys]).T[:-1]

def triangle(baseline):
    angles = np.linspace(0,2*np.pi,4)
    xs = 0.5773*baseline*np.sin(angles)
    ys = 0.5773*baseline*np.cos(angles)
    return np.array([xs,ys]).T[:-1]


def get_nuller_response():

    N = knull.make_nuller_mat5()

    telescope_array = pentagon2(baseline)

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

    k1 = response[1]-response[4]
    k2 = response[2]-response[3]

    return response, k1, k2 #return intensity per telescope

def get_nuller_response_tri():

    N = knull.make_nuller_mat3()

    telescope_array = triangle(baseline)

    fov = wavelength/baseline*scale_factor #?????
    sky_angles = np.linspace(-fov/2,fov/2,sz)

    xy = np.meshgrid(sky_angles, sky_angles, indexing='ij')

    #x,y are telescope positions in units of wavelength
    x = telescope_array[:,0]/wavelength
    y = telescope_array[:,1]/wavelength

    #Response is the 5 output electric fields as a function of the position on the sky
    response = np.zeros((3,sz,sz), dtype='complex')

    for i in range(3):
        for k in range(3):
            #Inputs have a phase equal to xy array - linear multiplied by spatial frequency
            #ul + vm?
            response[k] += np.exp(2*np.pi*1j*(xy[0]*x[i] + xy[1]*y[i]))*N[k,i] #

    response = np.abs(response)**2
    response /= (np.max(response[0])/3) #Normalise by flux per telescope

    #Create the kernel nulls. This is turning the output intensities into the kernel nulls (K in 2018 paper)

    k1 = response[1]-response[2]

    return response, k1#return intensity per telescope


r,k1,k2 = get_nuller_response()

rs = np.linspace(100,199,100)
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
plt.plot(rs*pix2mas/rad2mas*baseline/wavelength,k1_ave,label="k1")
plt.plot(rs*pix2mas/rad2mas*baseline/wavelength,k2_ave,label="k2")
#plt.axvline(x=HZAngle,color="k")
plt.legend()
plt.xlabel("Radial coordinate (mas)")
plt.ylabel("Transmission")
"""
r,k1 = get_nuller_response_tri()

rs = np.linspace(0.01,199)
k1_ave = []

for r in rs:
    k1_ave.append(azimuthal_rms(k1,r))

plt.figure(1)
plt.imshow(k1)
plt.figure(3)
plt.plot(rs*pix2mas,k1_ave,label="k1")
plt.axvline(x=HZAngle,color="k")
plt.legend()
plt.xlabel("Radial coordinate (mas)")
plt.ylabel("Transmission")
"""
