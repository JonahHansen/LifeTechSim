import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as const
from opticstools import knull

rad2mas = np.degrees(1)*3600e3 #Number of milliarcsec in one radian

wavelength = 14e-6 #m
scale_factor = 1.2
base_scale_factor = 1.04

sz = 400

#Tau Ceti
L = 0.2 #Lsol
Dist = 10.65 #Pc

HZAngle = np.sqrt(L)*1000/Dist

baseline = base_scale_factor*wavelength*rad2mas/HZAngle
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

from engine.nullers.four_telescopes import get_nuller_response


outputs = get_nuller_response(baseline,fov,sz,wavelength)

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
