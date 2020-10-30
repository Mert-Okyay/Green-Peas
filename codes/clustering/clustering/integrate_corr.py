import numpy as np
import scipy.integrate as intg
from scipy.stats import poisson
import time
from scipy.interpolate import InterpolatedUnivariateSpline as _spline
from scipy.integrate import simps


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
