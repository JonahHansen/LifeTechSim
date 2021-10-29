import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as const
from opticstools import knull

rad2mas = np.degrees(1)*3600e3 #Number of milliarcsec in one radian

"""
Calculate the RMS average over azimuthal angles for a given response map and
radial position

Inputs:
    image = map to average over
    r = radial position to average over angles

Outputs:
    RMS average
"""
def azimuthal_rms(image,r):
    n_angles = 1000
    angles = np.linspace(0,2*np.pi,n_angles)

    centre = (int(image.shape[0]/2),int(image.shape[1]/2))
    sum = 0
    for theta in angles:
        x = centre[0] + r*np.cos(theta)
        y = centre[1] + r*np.sin(theta)

        a = image[int(x),int(y)]

        sum += a**2

    return np.sqrt(sum/n_angles)

######################
architecture = 7

wavelength = 15e-6 #m
sz = 1500 #Size of grid
fov_scale_factor = 5 #Field of view scale factor

#What angle is the baseline to be optimised for?
L = 0.6 #Lsol
Dist = 6.65 #Pc
HZAngle = np.sqrt(L)*1000/Dist
HZAngle = 69.65

######################

#Set architecture, and define the baseline scale factor
if architecture == 1:
    from engine.nullers.bracewell import get_nuller_response
    architecture_verbose = "Bracewell four telescope nuller"
    base_scale_factor = 0.590

elif architecture == 2:
    from engine.nullers.linear import get_nuller_response
    architecture_verbose = "Linear four telescope nuller"

elif architecture == 3:
    from engine.nullers.three_telescopes import get_nuller_response
    architecture_verbose = "Three telescope kernel nuller"
    base_scale_factor = 0.666

elif architecture == 4:
    from engine.nullers.four_telescopes import get_nuller_response
    architecture_verbose = "Four telescope kernel nuller, optimised for K1"
    base_scale_factor = 0.4

elif architecture == 5:
    from engine.nullers.four_telescopes import get_nuller_response
    architecture_verbose = "Four telescope kernel nuller, optimised for K2"
    base_scale_factor = 0.4

elif architecture == 6:
    from engine.nullers.four_telescopes import get_nuller_response
    architecture_verbose = "Four telescope kernel nuller, optimised for K3"
    base_scale_factor = 0.4

elif architecture == 7:
    from engine.nullers.five_telescopes import get_nuller_response
    architecture_verbose = "Five telescope kernel nuller, optimised for adjacent telescopes (K1)"
    base_scale_factor = 1.028

elif architecture == 8:
    from engine.nullers.five_telescopes import get_nuller_response
    architecture_verbose = "Five telescope kernel nuller, optimised for diagonal telescopes (K2)"
    base_scale_factor = 0.660 #= approx 1.03*0.619 (where 0.619 is the conversion between a side and diagonal of a pentagon)


#baseline, fov defined as in the normal simulation
baseline = base_scale_factor*wavelength*rad2mas/HZAngle
fov = 2*fov_scale_factor*HZAngle/rad2mas
pix2mas = fov*rad2mas/sz

#Get response maps
outputs = get_nuller_response(baseline,fov,sz,wavelength)

rs = np.linspace(0.01,sz/2-0.01,100)

ks = []
i = 1
for (res,k) in outputs: #For each kernel output
    k_ave = []
    for r in rs:
        k_ave.append(azimuthal_rms(k,r)) #RMS azimuthal average for different radial positions
    ks.append(k_ave)
    plt.figure(i)
    plt.imshow(k) #Plot the kernel map
    i+=1

#Plot radial RMS average for each kernel output in units of lambda/B
plt.figure(i)
for j in range(len(outputs)):
    plt.plot(rs*pix2mas/rad2mas*baseline/wavelength,ks[j],label="K%s"%(j+1))
plt.title("Transmission per kernel output")
plt.ylabel("Transmission per telescope")
plt.xlabel(r"Angular distance ($\lambda/B$)")
plt.legend()

#Plot radial RMS average for the sum of each kernel output in units of lambda/B
ks = np.array(ks)
plt.figure(i+1)
plt.plot(rs*pix2mas/rad2mas*baseline/wavelength,np.sum(ks,axis=0))
plt.title("Total transmission")
plt.ylabel("Transmission per telescope")
plt.xlabel(r"Angular distance ($\lambda/B$)")

#Plot radial RMS average for each kernel output in units of mas
plt.figure(i+2)
for j in range(len(outputs)):
    plt.plot(rs*pix2mas,ks[j],label="K%s"%(j+1))
plt.title("Transmission per kernel output")
plt.ylabel("Transmission per telescope")
plt.xlabel(r"Angular distance (mas)")
