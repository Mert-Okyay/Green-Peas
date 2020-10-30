import numpy as np 
from numpy import random as rand
from scipy.stats import gaussian_kde
from astropy.table import Table, vstack
import matplotlib.pyplot as plt
import sys
from numpy.lib.recfunctions import append_fields
from astropy import units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord, match_coordinates_sky
from kde import weighted_gaussian_kde
from scipy import interpolate
from astropy.cosmology import FlatLambdaCDM,Planck15
import math
import time
from progressbar import Bar, AdaptiveETA, Percentage, ProgressBar
from clustering.utils import *

def simpleRand(data, bw=.2, dr=20, b=np.arange(1.4,3.06,0.01),cosmo=FlatLambdaCDM(H0=70, Om0=0.3)):
	'''
	keeps ra and dec within patch
	pulls from z distribution of origional patch
	'''
	nr=len(data)*dr
	#get smoothed z dist
	#get min and max for ra and dec of patch
	kernel = gaussian_kde(data['z'], bw_method=bw)
	evaluated=kernel.evaluate(b)

    #Create randoms
	z=[]
	r=[]
	d=[]
	q=np.cumsum(evaluated)/np.cumsum(evaluated)[-1]
	inv=interpolate.interp1d(q, b)

	ra_max = 150.645706
	ra_min = 149.710846
	dec_min = 1.718882
	dec_max = 2.63344
	pos = np.random.randint(len(data), size=nr)

	#progressbar
	widgets = [Percentage(),' ', Bar(),' ', AdaptiveETA()]
	pbar = ProgressBar(widgets=widgets, maxval=dr*len(data['z'])).start()
	for j in range(dr*len(data['z'])):
		a=rand.uniform(q.min(),q.max())
		z.append(inv(a))
		r.append(data[pos[j]]['ra'])
		d.append(data[pos[j]]['dec'])
		pbar.update(j+1)
	pbar.finish()
	temp=Table([r, d, z], names=('ra', 'dec', 'z'), meta={'name': 'first table'})
	random=np.array(temp)
	random=z_to_cdist(random,cosmo=cosmo)

	return random


