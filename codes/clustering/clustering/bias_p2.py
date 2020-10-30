from __future__ import division
import numpy as np
import math
from astropy import units as u
from astropy.coordinates import SkyCoord, match_coordinates_sky
import astropy.coordinates as coord
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table
from astropy.io import fits
import time
from numpy.lib.recfunctions import append_fields
from kde import weighted_gaussian_kde
from astropy.cosmology import FlatLambdaCDM,Planck15
from scipy.integrate import quad,romberg
import sys
from scipy.interpolate import interp1d

import matplotlib.pyplot as plt


from halomod import HaloModel, bias
from hmf import MassFunction

from scipy.interpolate import InterpolatedUnivariateSpline as _spline
from scipy.integrate import simps
from halomod.halo_model import HaloModel
from hmf._cache import cached_quantity, parameter
from halomod.halo_exclusion import dblsimps
from hmf.cosmo import Cosmology as csm
import warnings

print(wp_mm(sys.argv[1], sys.argv[2],Mmin=9,Mmax=16.5))

def wp_mm(z,pimax,Mmin=9,Mmax=16.5)
	hmf = MassFunction(Mmin=Mmin,Mmax=Mmax,z=0)
	wpdm = projected_corr_gal(hm.r, hm.corr_mm_lin, pimax, rp_out=None)
	return wpdm

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