import numpy as np 
from numpy.lib.recfunctions import append_fields
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import SkyCoord
from astropy import units as u
import sys
from scipy import interpolate

path = '/Users/justin/Astro/Yale18/Clustering/clustering/XMM_XXL/'
fluxarea=np.genfromtxt(path+'XMM_XXL_tot_af_0.5-10keV.txt', delimiter=',')
F = fluxarea[:,0]
A = fluxarea[:,1]
f = interpolate.interp1d(F,A)

def area(log_flux):
    if log_flux>5.990618755405745e-13:
        return 19.142640710580287999999
    elif log_flux<1.3595747888441975e-15:
        return 1e-6
    else:
        return f(log_flux)

def include_area_weights(cat):
    p_arr=[]
    Xarea=19.142640710580288
    for f in cat['flux']:
        p = area(f)/Xarea
        p_arr.append(1./p)
    norm = len(cat)/np.sum(p_arr)
    w_arr=norm*np.array(p_arr)
    if 'weight' in cat.dtype.names:
        w_arr = w_arr * cat['weight']
    temp = list(zip(cat['id'],cat['z'],cat['nh'],cat['lum'],cat['ra'],cat['dec'],cat['cdist'], cat['flux'],cat['hr'],w_arr))
    new = np.zeros((len(cat),), dtype=[('id', '<a8'),('z', '<f8'),('nh', '<i8'),('lum', '<f8'),('ra', '<f8'),('dec', '<f8'),('cdist', '<f8'),('flux', '<f8'),('hr', '<f8'),('weight', '<f8')])
    
    new[:] = temp
    return new

