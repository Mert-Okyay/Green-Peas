import numpy as np 
import random
from scipy.stats import gaussian_kde
from astropy.table import Table
import matplotlib.pyplot as plt
import sys
from numpy.lib.recfunctions import append_fields
from astropy import units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord, match_coordinates_sky
from kde import weighted_gaussian_kde
from scipy import interpolate

from clustering.utils import *

direc = '/Users/justin/Astro/Yale18/Clustering/clustering/stripe82x/'

def genrand(data,n,cosmo,width=.2,use_S82X_sens_map=True,plot=True,plot_filename=None):
	'''
	generates random catalog with random sky distribution and redshift
	To filter based on the S82X sensitivity map, set 'use_S82X_sens_map' to True
	
	'''
	goodz=data['z']>0
	if 'weight' in data.dtype.names:
		weights = data['weight'][goodz]
	else: weights=None
	z_arr = data['z'][goodz]
	ra_arr = data['ra'][goodz]
	dec_arr = data['dec'][goodz]
	ndata = len(np.unique(data['ra'][goodz]))
	
	#generate random redshifts
	if use_S82X_sens_map is True:
		n_rand = int(round(n*ndata*3))
	else: n_rand = int(round(n*ndata))
	z_grid = np.linspace(min(z_arr), max(z_arr), 1000)
	kde = weighted_gaussian_kde(z_arr, bw_method=width, weights=weights)
	kdepdfz=kde.evaluate(z_grid)
	zr_arr = generate_rand_from_pdf(pdf=kdepdfz, num=n_rand, x_grid=z_grid)
	
	#generate random sky coord in A013
	ra_max = 28.2
	ra_min = 14
	dec_min = -0.8
	dec_max = 0.8
	rar_arr = (np.random.rand(len(zr_arr))*(ra_max - ra_min))+ra_min
	#decrad = np.arcsin(np.random.rand(len(zr_arr))*2.-1) * u.radian
	#decr_arr = decrad.to(u.degree).value
	decr_arr = (np.random.rand(len(zr_arr))*(dec_max - dec_min))+dec_min
	
	
	temp = list(zip(zr_arr,rar_arr,decr_arr))
	rcat = np.zeros((len(zr_arr),), dtype=[('z', '<f8'),('ra', '<f8'),('dec', '<f8')])
	rcat[:] = temp
	
	if use_S82X_sens_map is True:
		if 'flux' not in data.dtype.names:
			print('no flux data in catalog found to filter based on sensitivity')
		else: 
			rcat = S82X_sensitivity_filter(data,rcat)

	randoms=rcat
	rcdists = np.array([cosmo.comoving_distance(z).value for z in randoms['z']])*cosmo.h
	randoms = append_fields(randoms, 'cdist', rcdists)
	randoms=np.array(randoms)
	
	print('number of randoms:', len(randoms))

	if plot:
		plot_zdist(data[goodz],randoms,z_grid,kdepdfz,plot_filename,weights=weights)

	return randoms

def S82X_sensitivity_filter(data,rcat):

	flux_arr = data['flux']
	n_rand = len(rcat)

	#get logN-logS distribution:
	fn = get_lognlogs()
	
	#generate random fluxes
	log_flux_grid = np.linspace(min(np.log10(flux_arr)), max(np.log10(flux_arr)), 1000)
	fpdf = fn(log_flux_grid)/np.sum(fn(log_flux_grid))
	log_rflux_arr = generate_rand_from_pdf(pdf=fpdf, num=n_rand, x_grid=log_flux_grid)
	fluxr_arr = 10**(log_rflux_arr)

	rcat = append_fields(rcat, 'flux', fluxr_arr)
	
	#filter based on sensitivity
	smap,wcs0 = get_S82Xsmap()
	good=[]
	for i,r in enumerate(rcat):
		ra=r['ra']
		dec=r['dec']
		flux = r['flux'] #ergs/s/cm^2
		py,px = wcs0.wcs_world2pix(ra,dec,1)
		px = np.round(px).astype(int)
		py = np.round(py).astype(int)
		sensitivity = smap[px,py]*10**-11
		#sensitivity map is a bit off for unknown reasons!! Fudge factor is used so that the data/random distributions match more closely
		if flux>.37*sensitivity:
			good = np.append(good,i)
	randoms=rcat[good.astype(int)]
	
	return randoms

def get_lognlogs():
	cat=np.genfromtxt(direc + 'lognlogs/s82_xmm_logn_logs0.5-10keV.txt')
	S = cat[:,0]
	N = cat[:,1]
	fN = interpolate.interp1d(np.log10(S),N)
	return fN

def get_S82Xsmap():
	'''Enter Galactic coordinates'''
	from astropy.wcs import WCS
	w0 = WCS(direc + 'exp_maps/ao13_bgmsk_expmap/Full/s82_0.5-10_sig5.1_corr.fits')
	h0 = fits.open(direc + 'exp_maps/ao13_bgmsk_expmap/Full/s82_0.5-10_sig5.1_corr.fits')
	s0 = h0[0].data
	smap = s0
	wcs = w0
	return smap,wcs

def generate_rand_from_pdf(pdf, num, x_grid):
	cdf = np.cumsum(pdf)
	cdf = cdf / cdf[-1]
	values = np.random.rand(num)
	value_bins = np.searchsorted(cdf, values)
	random_from_cdf = x_grid[value_bins]
	return random_from_cdf

def plot_zdist(data,randoms,z_grid,kdepdf,filename=None,weights=None):
	nd = len(data)
	nr = len(randoms)
	#if weights is None:
	#	plt.hist(data['z'], histtype='stepfilled', bins=20,alpha=0.2,color='k',density=True,label='data \n N='+str(nd))
	#else:
	plt.hist(data['z'], histtype='stepfilled', bins=20,alpha=0.2,color='k',density=True,label='data \n N='+str(nd),weights=weights)
	plt.hist(randoms['z'], histtype='step', bins=20,alpha=0.5,density=True,label='randoms\n N='+str(nr))
	plt.plot(z_grid, kdepdf, color='g', alpha=0.5, lw=3,label='smoothed PDF')
	plt.xlabel('z')
	plt.legend(loc='best', frameon=False)
	plt.tight_layout()
	if filename:
		plt.savefig(filename)