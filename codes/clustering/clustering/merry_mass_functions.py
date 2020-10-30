from __future__ import division
import numpy as np
import math
from astropy import units as u
from astropy.coordinates import SkyCoord, match_coordinates_sky
import astropy.coordinates as coord
from astropy.cosmology import FlatLambdaCDM,FLRW,Planck15
from astropy.table import Table
from astropy.io import fits
import time
from numpy.lib.recfunctions import append_fields
from kde import weighted_gaussian_kde
import sys
from scipy import interpolate, ndimage
from scipy.interpolate import interp1d
import scipy.integrate as intg
from scipy.stats import poisson
from scipy.interpolate import InterpolatedUnivariateSpline as _spline
from scipy.integrate import simps

from hmf import MassFunction

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

###for calculating wp(DM)
def power_to_corr_ogata(power, k, r, N=640, h=0.005):
    """
    Use Ogata's method for Hankel Transforms in 3D for nu=0 (nu=1/2 for 2D)
    to convert a given power spectrum to a correlation function.
    """
    lnk = np.log(k)
    spl = _spline(lnk, power)
    roots = np.arange(1, N + 1)
    t = h*roots
    s = np.pi*np.sinh(t)
    x = np.pi*roots*np.tanh(s/2)

    dpsi = 1 + np.cosh(s)
    dpsi[dpsi != 0] = (np.pi*t*np.cosh(t) + np.sinh(s))/dpsi[dpsi != 0]
    sumparts = np.pi*np.sin(x)*dpsi*x

    allparts = sumparts*spl(np.log(np.divide.outer(x, r))).T
    return np.sum(allparts, axis=-1)/(2*np.pi**2*r**3)

def projected_corr_gal(r, xir, pimax, rp_out=None):
    """
    Projected correlation function w(r_p).
    From Beutler 2011, eq 6.
    To integrate, we perform a substitution y = x - r_p.
    Parameters
    ----------
    r : float array
        Array of scales, in [Mpc/h]
    xir : float array
        Array of xi(r), unitless
    """

    if rp_out is None:
        rp_out = r

    lnr = np.log(r)
    lnxi = np.log(xir)

    p = np.zeros_like(rp_out)
    fit = _spline(r, xir, k=3)  # [self.corr_gal > 0] maybe?
    f_peak = 0.01
    a = 0

    for i, rp in enumerate(rp_out):
        #updated limit of integration (see Merry's paper)
        rlim = np.sqrt(pimax**2 + rp**2)

        if a != 1.3 and i < len(r) - 1:
            # Get log slope at rp
            ydiff = (lnxi[i + 1] - lnxi[i]) / (lnr[i + 1] - lnr[i])
            # if the slope is flatter than 1.3, it will converge faster, but to make sure, we cut at 1.3
            a = max(1.3, -ydiff)
            theta = _get_theta(a)

        min_y = theta * f_peak ** 2 * rp

        # Get the upper limit for this rp
        ylim = rlim - rp

        # Set the y vector for this rp
        y = np.logspace(np.log(min_y), np.log(ylim), 1000, base=np.e) 

        # Integrate
        integ_corr = fit(y + rp)
        integrand = (y + rp) * integ_corr / np.sqrt((y + 2 * rp) * y)
        p[i] = simps(integrand, y) * 2

    return p

def _get_theta(a):
    theta = 2 ** (1 + 2 * a) * (7 - 2 * a ** 3 + 3 * np.sqrt(5 - 8 * a + 4 * a ** 2) + a ** 2 * (9 + np.sqrt(5 - 8 * a + 4 * a ** 2)) -
                       a * (13 + 3 * np.sqrt(5 - 8 * a + 4 * a ** 2))) * ((1 + np.sqrt(5 - 8 * a + 4 * a ** 2)) / (a - 1)) ** (-2 * a)
    theta /= (a - 1) ** 2 * (-1 + 2 * a + np.sqrt(5 - 8 * a + 4 * a ** 2))
    return theta






    

# ## bias($\nu$) functions

#Tinker 2010
def bias_T10(nu):
    delta_halo = 200
    y = np.log10(delta_halo)
    delta_c = 1.686
    
    A = 1.0 + 0.24 * y * np.exp(-(4 / y) ** 4)
    a = 0.44 * y - 0.88
    C = 0.019 + 0.107 * y + 0.19 * np.exp(-(4 / y) ** 4)
    #nu = np.sqrt(self.nu)
    B = .183
    c = 2.4
    b = 1.5
    return 1 - A * nu ** a / (nu ** a + delta_c ** a) + B * nu ** b + C * nu ** c

#Tinker 2005
def bias_T05(nu):
    delta_c = 1.686
    a= .707
    b=.35
    c=.8
    #print(nu,delta(z))
    return 1+(1./(np.sqrt(a)*delta_c)) * ( (np.sqrt(a)*a*nu**2) + np.sqrt(a)*b*(a*nu**2)**(1-c) -                                            ((a*nu**2)**c/((a*nu**2)**c + b*(1-c)*(1-c/2.))) )

#Sheth 2001
def bias_S01(nu):
    delta_c = 1.686
    a= .707
    b=.5
    c=.6
    #print(nu,delta(z))
    return 1+(1./(np.sqrt(a)*delta_c)) * ( (np.sqrt(a)*a*nu**2) + np.sqrt(a)*b*(a*nu**2)**(1-c) -                                            ((a*nu**2)**c/((a*nu**2)**c + b*(1-c)*(1-c/2.))) )


# ## To obtain halo mass from bias:

# halo mass to nu
def m_to_nu(M,z,cosmo):
    delta_c = 1.686
    hmf = MassFunction(Mmin=9,Mmax=16,z=z,cosmo_model=Planck15)
    g = hmf.growth_factor
    #delta_c = 1.686/g
    return delta_c/sigma(M,z,cosmo)


#DM linear power spectrum, using transfer function from Eisenstein & Hu (1998) from hmf
def P_k(k,z,cosmo):
    hmf = MassFunction(Mmin=9,Mmax=16,z=z,cosmo_model=cosmo,transfer_model='EH_BAO')
    f = interp1d(hmf.k,hmf.power,kind='cubic')
    return f(k)

def mean_density0(cosmo):
    return (cosmo.Om0 * cosmo.critical_density0 / (cosmo.h**2)).to(u.Msun/u.Mpc**3).value

def radii(mass,cosmo):
    return (3.*mass / (4.*np.pi * mean_density0(cosmo))) ** (1. / 3.)

def k_space(kr):
    return np.where(kr>1.4e-6,(3 / kr ** 3) * (np.sin(kr) - kr * np.cos(kr)),1)

#root-mean square of mass density fluctuations withina sphere containing mass M
def sigma(mass,z,cosmo):
    lnk_min=np.log(1e-8)
    lnk_max=np.log(2e4)
    dlnk=0.05
    k_arr = np.exp(np.arange(lnk_min, lnk_max, dlnk))

    r = radii(mass,cosmo)
    rk = np.outer(r,k_arr)

    # we multiply by k because our steps are in logk.
    rest = P_k(k_arr,z,cosmo) * k_arr**3 
    integ = rest*k_space(rk)**2
    sig = (0.5/np.pi**2) * intg.simps(integ,dx=dlnk,axis=-1)
    return np.sqrt(sig)

def Mh(b,z,cosmo):
    m_arr = 10 ** np.arange(9, 16, 0.01)
    nu_arr = m_to_nu(m_arr,z,cosmo)
    f=interp1d(bias_T10(nu_arr),m_arr,kind='cubic')
    return f(b)

