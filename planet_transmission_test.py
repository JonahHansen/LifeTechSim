import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as const
from opticstools import knull

rad2mas = np.degrees(1)*3600e3 #Number of milliarcsec in one radian

wavelength = 14e-6 #m
scale_factor = 2
base_scale_factor = 1.03

sz = 400

#Tau Ceti
L = 0.6 #Lsol
Dist = 6.65 #Pc

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

from engine.nullers.five_telescopes import get_nuller_response

#1.69 for four telescopes
outputs = get_nuller_response(baseline,fov,sz,wavelength)

rs = np.linspace(0.01,199,100)

ks = []
i = 1
for (res,k) in outputs:
    k_ave = []
    for r in rs:
        k_ave.append(azimuthal_rms(k,r))
    ks.append(k_ave)
    plt.figure(i)
    plt.imshow(k)
    i+=1

plt.figure(i)
for j in range(len(outputs)):
    plt.plot(rs*pix2mas/rad2mas*baseline/wavelength,ks[j],label="K%s"%(j+1))
plt.title("Transmission per kernel output")
plt.ylabel("Transmission per telescope")
plt.xlabel(r"Angular distance ($\lambda/B$)")
plt.legend()


ks = np.array(ks)
plt.figure(i+1)
plt.plot(rs*pix2mas/rad2mas*baseline/wavelength,np.sum(ks,axis=0))
plt.title("Total transmission")
plt.ylabel("Transmission per telescope")
plt.xlabel(r"Angular distance ($\lambda/B$)")
