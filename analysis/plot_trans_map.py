sys.path.append("..")
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import cycler

rad2mas = np.degrees(1)*3600e3 #Number of milliarcsec in one radian

architecture = 10
base_wave = 15

arch_names = ["Bracewell","Kernel 3","Kernel 4","Kernel 5 (1.03)","Kernel 5 (0.66)","Kernel 5 (2.67)","Kernel 5 (1.68)"]

min_wave = 3
max_wave = 18
num_channels = 10
base_wave*=1e-6

sz = 2000 #Size of grid

#What angle is the baseline to be optimised for?
L = 0.6 #Lsol
Dist = 6.65 #Pc
HZAngle = np.sqrt(L)*1000/Dist
HZAngle = 32.77

######################
class Spectrograph():
    def __init__(self, wave_min, wave_max, baseline_wave, n_channels):

        #Bandwidth of spectrograph
        self.bandwidth = (wave_max-wave_min)*1e-6
        #Mean wavelength of spectrograph
        self.mean = (wave_max+wave_min)/2*1e-6
        #Min and max wavelengths
        self.wave_min = wave_min*1e-6
        self.wave_max = wave_max*1e-6
        #Number of channels
        self.n_channels = n_channels

        #Wavelength to set the baseline (optimisation baseline)
        self.baseline_wave = baseline_wave*1e-6

        #Channel borders (start of each channel)
        self.channel_borders = np.linspace(self.wave_min,self.wave_max,n_channels+1)[:-1]
        #Size of channel
        self.dlambda = self.channel_borders[1]-self.channel_borders[0]
        #Channel centres (middle of each channel)
        self.channel_centres = (self.channel_borders + self.dlambda/2)
        #Effective resolution of the spectrograph
        self.eff_resolution = self.mean/self.dlambda

spec = Spectrograph(min_wave,max_wave,base_wave,num_channels)
fov_scale_factor = base_wave/(spec.channel_centres[0]) + 0.1

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
    architecture_verbose = "Four telescope kernel nuller"
    base_scale_factor = 0.4

elif architecture == 7:
    from engine.nullers.five_telescopes import get_nuller_response
    architecture_verbose = "Five telescope kernel nuller, optimised at 1.03 B/lambda"
    base_scale_factor = 1.028

elif architecture == 8:
    from engine.nullers.five_telescopes import get_nuller_response
    architecture_verbose = "Five telescope kernel nuller, optimised at 0.66 B/lambda"
    base_scale_factor = 0.660 #= approx 1.03*0.619 (where 0.619 is the conversion between a side and diagonal of a pentagon)

elif architecture == 9: #THESE ARE ALTERNATIVES - MAY BE BETTER!
    from engine.nullers.five_telescopes import get_nuller_response
    architecture_verbose = "Five telescope kernel nuller, optimised at 2.67 B/lambda"
    base_scale_factor = 2.67

elif architecture == 10:
    from engine.nullers.five_telescopes import get_nuller_response
    architecture_verbose = "Five telescope kernel nuller, optimised at 1.68 B/lambda"
    base_scale_factor = 1.68 #= approx 1.03*0.619 (where 0.619 is the conversion between a side and diagonal of a pentagon)

def azimuthal_array(image,r):

    n_angles = 5000
    angles = np.linspace(0,2*np.pi,n_angles)

    centre = (int(image.shape[0]/2),int(image.shape[1]/2))

    #Planet out of field of view!
    if r > image.shape[0]/2:
        return 0

    arr = []
    pos = []
    for theta in angles:
        x = centre[0] + r*np.cos(theta)
        y = centre[1] + r*np.sin(theta)

        a = image[int(x),int(y)]

        arr.append(a)
        pos.append(np.array([x,y]))

    return np.array(arr), pos

def round_sig_figs(x, p):
    x_positive = np.where(np.isfinite(x) & (x != 0), np.abs(x), 10**(p-1))
    mags = 10 ** (p - 1 - np.floor(np.log10(x_positive)))
    return np.round(x * mags) / mags

#baseline, fov defined as in the normal simulation
baseline = base_scale_factor*base_wave*rad2mas/HZAngle
fov = 2*fov_scale_factor*HZAngle/rad2mas
pix2mas = fov*rad2mas/sz
wave_pix2mas = pix2mas/base_wave

#Get response maps
outputs = get_nuller_response(baseline,fov,sz,base_wave)

signal = []
pos = []
for (res,k) in outputs: #For each kernel output
    temp_signal = []
    temp_pos = []
    for wave in spec.channel_centres:
        planet_pos = HZAngle/(wave_pix2mas*wave) #in pixels

        #Characterisation mode is maximum of kernel azimuth
        #import pdb; pdb.set_trace()
        p_trans,p_pos = azimuthal_array(k,planet_pos)
        temp_signal.append(np.abs(p_trans))
        temp_pos.append(p_pos)

    signal.append(np.array(temp_signal))
    pos.append(np.array(temp_pos))


signal = np.array(signal)
pos = np.array(pos)
summed_signal = np.sum(signal,axis=(0,1))
arg = np.argmax(summed_signal)
pos2 = pos[:,:,arg]


color = plt.cm.winter(np.linspace(0, 1,num_channels))
mpl.rcParams['axes.prop_cycle'] = cycler.cycler('color', color)

i=0
for (res,k) in outputs:
    plt.figure(i)
    plt.clf()
    plt.imshow(k,cmap="RdGy",extent=[-fov_scale_factor,fov_scale_factor,-fov_scale_factor,fov_scale_factor]) #Plot the kernel map
    plt.plot(0,0,marker="*",c="orange",markersize=20)
    plt.colorbar(label="Transmission per telescope flux")
    pos_k = pos2[i]
    for j in range(len(pos_k)):
        plt.scatter((pos_k.T[1][j]-sz/2)/(sz/2)*fov_scale_factor,(pos_k.T[0][j]-sz/2)/(sz/2)*fov_scale_factor,label=round_sig_figs(spec.channel_centres[j]*1e6,5))
    i+=1
    plt.legend()
    plt.xlabel("Angular position normalised by\n baseline optimised position")
    plt.ylabel("Angular position normalised by\n baseline optimised position")
    plt.title("K%s Transmission map and planet positions as a function of wavelength \n for the %s architecture\n and a reference wavelength of %s um" %(i,architecture_verbose,round_sig_figs(base_wave*1e6,2)))

plt.show()
