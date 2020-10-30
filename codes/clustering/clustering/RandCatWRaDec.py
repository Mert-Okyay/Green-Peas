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

def randBlockCat(data, bw=.2, dr=20, b=np.arange(0.1,7.01,0.01),cosmo=FlatLambdaCDM(H0=70, Om0=0.3)):
	'''
	keeps ra and dec within patch
	pulls from z distribution of origional patch
	'''

	#get smoothed z dist
	#get min and max for ra and dec of patch
	kernel = gaussian_kde(data['z'], bw_method=bw)
	evaluated=kernel.evaluate(b)

    #Create randoms
	r=[]
	d=[]
	z=[]
	q=np.cumsum(evaluated)/np.cumsum(evaluated)[-1]
	inv=interpolate.interp1d(q, b)

	#progressbar
	widgets = [Percentage(),' ', Bar(),' ', AdaptiveETA()]
	pbar = ProgressBar(widgets=widgets)
	for j in pbar(range(dr*len(data['z']))):
		rinp=rand.randint(0,len(data['z'])-1)
		r.append(data['ra'][rinp])
		d.append(data['dec'][rinp])
		a=rand.uniform(q.min(),q.max())
		z.append(inv(a))
		pbar.update(j+1)
	temp=Table([r, d, z], names=('ra', 'dec', 'z'), meta={'name': 'first table'})
	temp[temp['ra']<0]['ra']+=360
	random=np.array(temp)
	random=z_to_cdist(random,cosmo=cosmo)

	return random

def randPointBlockCat(patches, bw=.2, dr=20, b=np.arange(0.1,7.01,0.01),cosmo=FlatLambdaCDM(H0=70, Om0=0.3)):
	'''
	keeps ra and dec within patch
	pulls from z distribution of origional patch
	'''

	#get smoothed z dist of patches
	#get min and max for ra and dec of patch
	evalu=[]
	dmm=[]
	rmm=[]
	for p in patches:
		if (len(p['z'])>1):
			kernel = gaussian_kde(p['z'], bw_method=bw)
			evaluated=kernel.evaluate(b)
			evalu.append(evaluated)
			dmm.append([p['dec'].min(),p['dec'].max()])
			rmm.append([p['ra'].min(),p['ra'].max()])
		else:
			evalu.append(None)

    #Create randoms for each patch
	random=[]
	for i, p in enumerate(patches):
		if(len(p['z'])>1):
		    r=[]
		    d=[]
		    z=[]
		    q=np.cumsum(evalu[i])/np.cumsum(evalu[i])[-1]
		    inv=interpolate.interp1d(q, b)

		    #progressbar
		    widgets = [Percentage(),' ', Bar(),' ', AdaptiveETA()]
		    pbar = ProgressBar(widgets=widgets)

		    for j in pbar(range(dr*len(p['z']))):
		        r.append(rand.uniform(rmm[i][0],rmm[i][1]))
		        d.append(rand.uniform(dmm[i][0],dmm[i][1]))
		        a=rand.uniform(q.min(),q.max())
		        z.append(inv(a))
		        pbar.update(j+1)
		    temp=Table([r, d, z], names=('ra', 'dec', 'z'), meta={'name': 'first table'})
		    temp[temp['ra']<0]['ra']+=360
		    random.append(temp)

    #combine randomized patches for a total random catalog
	total=Table([[],[],[]], names=('ra', 'dec', 'z'), meta={'name': 'first table'})
	for i, r in enumerate(random):
	    total=vstack([total,r])

	temp=total[total['ra']<0]
	for t in temp:
		t['ra']+=360
	total[total['ra']<0]=temp
	total=np.array(total)
	total=z_to_cdist(total,cosmo=cosmo)

	return total


def XMMXXLrand(data, bw=.2, dr=20, b=np.arange(0.1,7.01,0.01),cosmo=FlatLambdaCDM(H0=70, Om0=0.3)):
	'''
	Generates random catalog by splitting into patches based on programname
	randomizes ra and dec within patch
	pulls from z distribution of origional patch
	'''
	patches=[]
	patches.append(data[data['programname']==b'boss'])
	patches.append(data[(data['programname']==b'eboss')&(data['ra']>34)&(data['ra']<38.25)&(data['dec']>-6)&(data['dec']<-3.75)])
	patches=np.array(patches)
	return randPointBlockCat(patches,bw,dr,b,cosmo)

def XMMagnRand(agn, bw=.2, dr=20, b=np.arange(0.1,7.01,0.01),cosmo=FlatLambdaCDM(H0=70, Om0=0.3)):
	'''
	Generates random catalog by splitting into patches based on programname
	randomizes ra and dec within patch
	pulls from z distribution of origional patch
	'''
	patches=[]
	patches.append(agn[agn['dec']<-7])
	patches.append(agn[(agn['dec']>=-7)&(agn['dec']<-6.5)])
	patches.append(agn[(agn['dec']>=-6.5)&(agn['dec']<-5.92)])
	patches.append(agn[(agn['dec']>=-5.92)&(agn['dec']<-5.4)])
	patches.append(agn[(agn['dec']>=-5.4)&(agn['dec']<-4.9)])
	patches.append(agn[(agn['dec']>=-4.9)&(agn['dec']<-4.688)])
	patches.append(agn[(agn['dec']>=-4.688)&(agn['dec']<-3.95)])
	patches.append(agn[(agn['dec']>=-3.95)&(agn['dec']<-3.25)])
	patches.append(agn[(agn['dec']>=-3.25)])
	patches=np.array(patches)
	return randPointBlockCat(patches,bw,dr,b,cosmo)



