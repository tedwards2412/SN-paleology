import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import swordfish as sf
from random import *
from tqdm import tqdm
import scipy.interpolate as interpolate
from scipy.optimize import minimize
from scipy.integrate import quad
import paleopy as paleopy
from WIMpy import DMUtils as DMU
from scipy.special import gamma
from scipy.special import erf

import os
R_E = 8.12 # kpc - Table A1 in 1807.09409 
# R_E = 8.7 # kpc - FIXME: Should this be changed?


f_Earth = None #Heliocentric PDF for SN distances

# Fermi-Dirac spectrum
def dNdE(E, Etot, T, alpha):
    numer = Etot*((1+alpha)**(1+alpha))*(E**alpha)*np.exp(-(1+alpha)*E/T)
    return numer/gamma(1+alpha)/(T**(2+alpha))

def inverse_transform_sampling(function, x_range, nbins=100, n_samples=1000):
    bins = np.linspace(x_range[0], x_range[-1], num=nbins)
    pdf = function(np.delete(bins,-1) + np.diff(bins)/2)
    Norm = np.sum(pdf*np.diff(bins))
    pdf /= Norm
    cum_values = np.zeros(bins.shape)
    cum_values[1:] = np.cumsum(pdf*np.diff(bins))
    inv_cdf = interpolate.interp1d(cum_values, bins)
    r = np.random.rand(n_samples)
    return inv_cdf(r)

def calc_f():
    global f_Earth
    
    if (f_Earth is None):
        
        R_gal = 50
        N_SN = int(1e5)

        def calc_dist(x, y, z):
            x_E = R_E
            y_E = 0.0
            z_E = 0.0
            return np.sqrt((x-x_E)**2+(y-y_E)**2+(z-z_E)**2)

        phi = np.random.uniform(0,2*np.pi, size=N_SN)
        z = np.random.exponential(scale=0.33, size=N_SN)

        l_c = 2.5 #kpc
        r_B = 2.0 #kpc
        sigma_0 = 611e6 #Msun/kpc^2

        def sigma_disc(r):
            #return r*sigma_0*l_c*((r-r_B)**2 + l_c**2)**-0.5
            r_d = 2.9 #kpc
            return r*np.exp(-r/r_d)

        x_temp = np.linspace(0,R_gal,num=100)
        r = inverse_transform_sampling(sigma_disc, x_temp, nbins=100, n_samples=N_SN)

        theta = np.pi/2.
        x = r * np.sin(theta) * np.cos(phi)
        y = r * np.sin(theta) * np.sin(phi)

        for i in tqdm(range(N_SN)):
            R = calc_dist(x, y, z)

        b_edge = np.linspace(0,R_gal,100)
        b_c = b_edge[:-1] + np.diff(b_edge)/2
        dist_bin = np.histogram(R, bins=b_edge, normed=True)

        f_Earth = interpolate.interp1d(b_c, dist_bin[0], bounds_error=False, fill_value=(0.0,0.0))
        
        xlist = np.linspace(0, 50, 100)
        plt.figure()
        plt.plot(xlist, f_Earth(xlist))
        plt.xlabel(r"Helio-centric distance $D$ [kpc]")
        plt.ylabel(r"$P(D)$ [kpc$^{-1}$]")
        plt.show()
    
    return f_Earth

def galactic_SN(E):
    T_v = np.array([13.3e3,14.6e3,15.0e3,15.0e3,15.0e3,15.0e3]) # keV
    Etot = np.array([6.0e52,4.3e52,2.0e52,2.0e52,2.0e52,2.0e52]) # ergs
    alpha = np.array([3., 3.3, 3., 3., 3., 3.])
    Etot *= 6.242e+8 # ergs to keV
    F = calc_f()

    
    dist_weight = lambda R: F(R)/4/np.pi/((R)**2.) 
    dndE = np.zeros_like(E)
    for i in range(6):
        dndE += dNdE(E, Etot[i], T_v[i], alpha[i])
    dNdE_gal = dndE*quad(dist_weight, 0., 30)[0]/(3.086e+21**2) # Convert kpc^-2 to cm^-2
    # returns in s^-1 keV^-1 cm^-2 - assuming one CC SN per second
    return dNdE_gal
    
def window(x, x_a, x_b, sigma):
    return 0.5*(erf((x - x_a)/(np.sqrt(2.0)*sigma)) - erf((x - x_b)/(np.sqrt(2.0)*sigma)))

def calcdRdx():   
    Epso = paleopy.Mineral("Epsomite")
    Epso.showProperties()
    Epso.showSRIM()

    Elist = np.logspace(-1,6,1000) # keV
    dndE = galactic_SN(Elist)

    dndE_gal = interpolate.interp1d(Elist*1e-3, dndE*1e3, bounds_error=False, fill_value=(0.0,0.0)) # converting keV --> MeV

    x0 = 15.0/2.0
    x = np.logspace(0,3,num=200)

    dRdx_temp = np.zeros_like(x)
    for i, nuc in enumerate(Epso.nuclei):
        if (nuc != "H"):
            xtemp = Epso.Etox_nuclei[nuc](Elist)
            dRdx_nuc = (np.vectorize(DMU.dRdE_CEvNS)(Elist, Epso.N_p[i], Epso.N_n[i], flux_name="user", flux_func = dndE_gal)
                                                *Epso.dEdx_nuclei[nuc](xtemp))
            temp_interp = interpolate.interp1d(xtemp, dRdx_nuc, fill_value='extrapolate')
            dRdx_temp += Epso.ratio_nuclei[nuc]*temp_interp(x)

    dRdx = dRdx_temp*1e6*365+1e-20 # kg/Myr/nm
    dRdx_interp = interpolate.interp1d(x, dRdx, fill_value='extrapolate')

    sigma = 15.
    x = np.logspace(np.log10(x0), 3, 70)

    x_c = x[:-1] + np.diff(x)/2
    dRdx_binned = np.zeros_like(x_c)

    for i in tqdm(range(x_c.size)):
        x1 = x_c[i] - 5.0*sigma
        x2 = x_c[i] + 5.0*sigma
        x1 = np.clip(x1, 0.1, 1e5)
        intge = lambda y: dRdx_interp(y)*window(y, x[i], x[i+1], sigma)
        dRdx_binned[i] = quad(intge, x1, x2)[0] + 1e-30

    return dRdx_binned

def calc_signal(mineral):
    Elist = np.logspace(-1,6,1000) # keV
    dndE = galactic_SN(Elist)
    dndE_gal = interpolate.interp1d(Elist*1e-3, dndE*1e3, bounds_error=False, fill_value=(0.0,0.0)) 

    x0 = 15.0/2.0
    x = np.logspace(0,3,num=200)

    dRdx_temp = np.zeros_like(x)
    for i, nuc in enumerate(mineral.nuclei):
        if (nuc != "H"):
            xtemp = mineral.Etox_nuclei[nuc](Elist)
            dRdx_nuc = (np.vectorize(DMU.dRdE_CEvNS)(Elist, mineral.N_p[i], mineral.N_n[i], flux_name="user", flux_func = dndE_gal)
                                                *mineral.dEdx_nuclei[nuc](xtemp))
            temp_interp = interpolate.interp1d(xtemp, dRdx_nuc, fill_value='extrapolate')
            dRdx_temp += mineral.ratio_nuclei[nuc]*temp_interp(x)

    dRdx = dRdx_temp*1e6*365+1e-20 # kg/Myr/nm
    dRdx_interp = interpolate.interp1d(x, dRdx, fill_value='extrapolate')
    
    sigma = 15.
    x = np.logspace(np.log10(x0), 3, 70)
    x_c = x[:-1] + np.diff(x)/2
    dRdx_binned = np.zeros_like(x_c)

    for i in range(x_c.size):
        x1 = x_c[i] - 5.0*sigma
        x2 = x_c[i] + 5.0*sigma
        x1 = np.clip(x1, 0.1, 1e5)
        intge = lambda y: dRdx_interp(y)*window(y, x[i], x[i+1], sigma)
        dRdx_binned[i] = quad(intge, x1, x2)[0] + 1e-30
    
    return dRdx_binned