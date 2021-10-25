#Old function graveyard



#Calculate stellar leakage flux through Mike's method
def Mike_stellar_leakage(star,response,pix2mas):
    #Find the null function numerically
    #Integral over theta, divide by two pi, as a function of radius
    r_ix,y2 = azimuthalAverage(response[1], returnradii=True, binsize=0.8)
    r_ix,y4 = azimuthalAverage(response[2], returnradii=True, binsize=0.8)
    #averaging response function over theta as a function of radius. yi is average value

    #Dodgy fit for a quartic and parabola
    second_order_coeff = np.median(y2[1:16]/r_ix[1:16]**2)/pix2mas**2
    fourth_order_coeff = np.median(y4[1:16]/r_ix[1:16]**4)/pix2mas**4

    r = np.linspace(0,1,200)
    I = sf.limb_darkening(r)

    #Total leakage is an integral of coeff multiplied by r^2 or r^4 and limb darkening intensity
    mn_r2 = (np.trapz(I*r**3, r)/np.trapz(I*r, r))**.5
    mn_r4 = (np.trapz(I*r**5, r)/np.trapz(I*r, r))**.25

    #convert to flux...
    leakage_2nd = second_order_coeff*(mn_r2*star.angRad)**2*star.flux
    leakage_4th = fourth_order_coeff*(mn_r4*star.angRad)**4*star.flux

    return leakage_2nd, leakage_4th


#How many multiples do I need of the baseline to keep the baseline between 10 and 600m
def baseline_checker(baseline):
    if baseline > 600:
        return 600, 0
    elif baseline >= 10:
        return baseline, 1
    else:
        n = 1
        while baseline < 10:
            n += 1
            baseline *= (2*n-1)
        return baseline, n
