from __future__ import division
import numpy as np
import math
from astropy import units as u
from astropy.coordinates import SkyCoord, match_coordinates_sky
import astropy.coordinates as coord
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table, join, vstack
from astropy.io import fits
import time
from numpy.lib.recfunctions import append_fields
from astropy.cosmology import FlatLambdaCDM,Planck15
import sys
from scipy import interpolate, stats
from scipy.integrate import quad as Integrate
from hmf import MassFunction

import scipy.integrate as intg
from scipy.stats import poisson
from scipy.interpolate import InterpolatedUnivariateSpline as _spline
from scipy.integrate import simps

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
import clustering.projected_correlation_functions as w
from clustering.zCosmosDEEP_randoms import simpleRand as rand
from clustering.zCosmosDEEP_utils import *
from clustering.utils import z_to_cdist
from clustering import merry_mass_functions as mmf

import clustering.halo_mass as hm

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import matplotlib
cmap = matplotlib.cm.get_cmap('viridis')


def effective(redshifts):
    total=0
    count=0
    for x, i in enumerate(redshifts):
        for y, j in enumerate(redshifts[x+1:]):
            total = total + (i+j)/2
            count = count + 1
    return total/count 

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

def projected_corr_gal(r, xir, rlim, rp_out=None):
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

def pipe(data, randoms, names, gal, gr, CF='X', pimax=60, nbins=10, m=25):
    n=len(data)
    cvalues=np.arange(0,1+1/(n+1),1/(n+1))
    c= []
    for v in cvalues:
        c.append(cmap(v))
    bins = np.logspace(-1, np.log10(pimax), nbins + 1)
    rp=[]
    wp=[]
    wp_err=[]
    cov=[]
    wpC=[]
    wpC_err=[]
    zeff=[]
    rpA, wpA, wpA_err, covA = w.auto_wp(data=gal,randoms=gr, bins=bins, pimax=pimax, m=m, estimator='L',survey='zCosDEEP')
    for i in range(n):
        rpX, wpX, wpX_err, covX = w.cross_wp(d1=data[i], r1=randoms[i], d2=gal,r2=gr, bins=bins, pimax=pimax, m=m, estimator='L',survey='zCosDEEP')
        rp.append(rpX)
        wp.append(wpX)
        wp_err.append(wpX_err)
        cov.append(covX)
        wpC.append(wp[i]*wp[i]/wpA)
        wpC_err.append(np.sqrt(((wp_err[i]/wp[i])**2)+((wp_err[i]/wp[i])**2)+((wpA_err/wpA)**2))*wpC[i])
        zeff.append(effective(data[i]['z']))
    hmf=MassFunction(n=1, z=zeff[0], Mmin=agn['logM'].min(), Mmax=agn['logM'].max())
    r= np.logspace(-1, 2, 1000)
    xir=power_to_corr_ogata(hmf.power,hmf.k,r)
    wpDM=projected_corr_gal(r, xir, pimax)

    if ((CF=='X')|(CF=='B')):
        XCF1,(ax1,ax2) = plt.subplots(1,2,figsize=(12,4))
        ax1.set_yscale("log", nonposy='clip')
        ax1.set_xscale("log", nonposx='clip')
        ax1.plot(r, wpDM, color=c[0], label='Dark Matter')
        for i in range(n):
            ax1.errorbar(rp[i],wp[i],yerr=wp_err[i],fmt='o',color=c[i+1], label=names[i]+' XCF')
            ax2.plot(rp[i],wp_err[i]/wp[i],color=c[i+1])
        ax1.errorbar(rpA,wpA,yerr=wpA_err,fmt='o',color=c[-1], label='Galaxies ACF')
        ax2.plot(rpA,wpA_err/wpA,color=c[-1])
        ax1.set_xlabel('r$_{p}$',fontsize=16) 
        ax1.set_ylabel('w$_{p}$',fontsize=16)
        ax1.legend(frameon=False,fontsize=11)
        ax1.set_title('zCosmos-DEEP')
        ax1.set_ylim(1,1000)
        ax2.set_yscale("log", nonposy='clip')
        ax2.set_xscale("log", nonposx='clip')
        ax2.set_xlabel('r$_{p}$',fontsize=16)
        ax2.set_ylabel('$\sigma_{p}/w_{p}$',fontsize=16)
        XCF1.show()
        XCF1.savefig('XCF1.pdf')

        XCF2,(ax1,ax2) = plt.subplots(1,2,figsize=(12,4))
        ax1.set_yscale("log", nonposy='clip')
        ax1.set_xscale("log", nonposx='clip')
        ax1.plot(r, wpDM, color=c[0], label='Dark Matter')
        for i in range(n):
            ax1.errorbar(rp[i],wp[i],yerr=wp_err[i],fmt='o',color=c[i+1], label=names[i]+' XCF')
            ax2.plot(rp[i],wp_err[i]/wp[i],color=c[i+1])
        ax1.errorbar(rpA,wpA,yerr=wpA_err,fmt='o',color=c[-1], label='Galaxies ACF')
        ax2.plot(rpA,wpA_err/wpA,color=c[-1])
        ax1.set_xlabel('r$_{p}$',fontsize=16)
        ax1.set_ylabel('w$_{p}$',fontsize=16)
        ax1.legend(frameon=False,fontsize=11)
        ax1.set_title('zCosmos-DEEP')
        ax1.set_ylim(1,1000)
        ax1.axis([0.8,15,1,1000])
        ax2.set_yscale("log", nonposy='clip')
        ax2.set_xscale("log", nonposx='clip')
        ax2.set_xlabel('r$_{p}$',fontsize=16)
        ax2.set_ylabel('$\sigma_{p}/w_{p}$',fontsize=16)
        ax2.axis([0.8,15,.05,4])
        XCF2.show()
        XCF2.savefig('XCF2.pdf')

    if ((CF=='A')|(CF=='B')):
        ACF1,(ax1,ax2) = plt.subplots(1,2,figsize=(12,4))
        ax1.set_yscale("log", nonposy='clip')
        ax1.set_xscale("log", nonposx='clip')
        ax1.plot(r, wpDM, color=c[0], label='Dark Matter')
        for i in range(n):
            ax1.errorbar(rpA,wpC[i],yerr=wpC_err[i],fmt='o',color=c[i+1], label=names[i]+' ACF')
            ax2.plot(rpA,wpC_err[i]/wpC[i],color=c[i+1])
        ax1.errorbar(rpA,wpA,yerr=wpA_err,fmt='o',color=c[-1], label='Galaxies ACF')
        ax2.plot(rpA,wpA_err/wpA,color=c[-1])
        ax1.set_xlabel('r$_{p}$',fontsize=16) 
        ax1.set_ylabel('w$_{p}$',fontsize=16)
        ax1.legend(frameon=False,fontsize=11)
        ax1.set_title('zCosmos-DEEP')
        ax1.set_ylim(1,1000)
        ax2.set_yscale("log", nonposy='clip')
        ax2.set_xscale("log", nonposx='clip')
        ax2.set_xlabel('r$_{p}$',fontsize=16)
        ax2.set_ylabel('$\sigma_{p}/w_{p}$',fontsize=16)
        ACF1.show()
        ACF1.savefig('ACF1.pdf')

        ACF2,(ax1,ax2) = plt.subplots(1,2,figsize=(12,4))
        ax1.set_yscale("log", nonposy='clip')
        ax1.set_xscale("log", nonposx='clip')
        ax1.plot(r, wpDM, color=c[0], label='Dark Matter')
        for i in range(n):
            ax1.errorbar(rpA,wpC[i],yerr=wpC_err[i],fmt='o',color=c[i+1], label=names[i]+' ACF')
            ax2.plot(rpA,wpC_err[i]/wpC[i],color=c[i+1])
        ax1.errorbar(rpA,wpA,yerr=wpA_err,fmt='o',color=c[-1], label='Galaxies ACF')
        ax2.plot(rpA,wpA_err/wpA,color=c[-1])
        ax1.set_xlabel('r$_{p}$',fontsize=16)
        ax1.set_ylabel('w$_{p}$',fontsize=16)
        ax1.legend(frameon=False,fontsize=11)
        ax1.set_title('zCosmos-DEEP')
        ax1.set_ylim(1,1000)
        ax1.axis([0.8,15,1,1000])
        ax2.set_yscale("log", nonposy='clip')
        ax2.set_xscale("log", nonposx='clip')
        ax2.set_xlabel('r$_{p}$',fontsize=16)
        ax2.set_ylabel('$\sigma_{p}/w_{p}$',fontsize=16)
        ax2.axis([0.8,15,.05,4])
        ACF2.show()
        ACF2.savefig('ACF2.pdf')










        