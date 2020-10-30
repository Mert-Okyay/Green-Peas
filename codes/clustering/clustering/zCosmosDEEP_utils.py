import numpy as np 
from numpy.lib.recfunctions import append_fields
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import SkyCoord
from astropy import units as u
import sys
from scipy import interpolate

path = '/Users/justin/Box/Astro/Yale18/Clustering/clustering/stripe82x/'
fluxarea=np.genfromtxt(path+'s82x_tot_af_0.5-10keV.txt')
F = fluxarea[:,0]
A = fluxarea[:,1]
f = interpolate.interp1d(F,A)

def area(log_flux):
    if log_flux>-13:
        return 31.288699999999999
    elif log_flux<-14.78:
        return 1e-6
    else:
        return f(log_flux)

def include_area_weights(cat):
    p_arr=[]
    s82area=31.2887
    for f in cat['flux']:
        p = area(np.log10(f))/s82area
        p_arr.append(1./p)
    norm = len(cat)/np.sum(p_arr)
    w_arr=norm*np.array(p_arr)
    if 'weight' in cat.dtype.names:
        w_arr = w_arr * cat['weight']
    if (('id' in cat.dtype.names)&('hr' in cat.dtype.names)&('ztype' in cat.dtype.names)&('nh' in cat.dtype.names)):
        temp = list(zip(cat['id'],cat['ztype'],cat['z'],cat['ra'],cat['dec'],cat['flux'],cat['hr'],cat['nh'],w_arr))
        new = np.zeros((len(cat),), dtype=[('id', '<f8'),('ztype', '<i8'),('z', '<f8'),('ra', '<f8'),('dec', '<f8'),('flux', '<f8'),('hr', '<f8'),('nh', '<f8'),('weight', '<f8')])
    elif (('id' in cat.dtype.names)&('hr' in cat.dtype.names)&('ztype' in cat.dtype.names)):
        temp = list(zip(cat['id'],cat['ztype'],cat['z'],cat['ra'],cat['dec'],cat['flux'],cat['hr'],w_arr))
        new = np.zeros((len(cat),), dtype=[('id', '<f8'),('ztype', '<i8'),('z', '<f8'),('ra', '<f8'),('dec', '<f8'),('flux', '<f8'),('hr', '<f8'),('weight', '<f8')])
    else:
        temp = list(zip(cat['z'],cat['ra'],cat['dec'],cat['flux'],w_arr))
        new = np.zeros((len(cat),), dtype=[('z', '<f8'),('ra', '<f8'),('dec', '<f8'),('flux', '<f8'),('weight', '<f8')])
    new[:] = temp
    return new