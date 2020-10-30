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

from clustering.utils import *

def genProgRand(data, bw=.2, dr=20, b=np.arange(0.1,2.01,0.01),cosmo=FlatLambdaCDM(H0=70, Om0=0.3)):
	'''
	Generates random catalog by splitting into patches based on programname
	randomizes ra and dec within patch
	pulls from z distribution of origional patch
	'''
	
	temp=data[data['ra']>180]
	for t in temp:
		t['ra']-=360
	data[data['ra']>180]=temp

	#define patches
	patches=np.empty(7, dtype=object)
	ebossALL=data[(data['programname']==b'eboss')]
	patches[0]=ebossALL[(ebossALL['ra']<=-3)]
	patches[1]=ebossALL[((ebossALL['ra']>-3)&(ebossALL['ra']<=0)&(ebossALL['dec']>-1)&(ebossALL['dec']<=1))]
	patches[2]=ebossALL[(((ebossALL['ra']>0)&(ebossALL['ra']<=14))|((((ebossALL['ra']>=-3)&(ebossALL['ra']<=0)))&((ebossALL['dec']>1)|(ebossALL['dec']<-1))))]
	patches[3]=ebossALL[((ebossALL['ra']>14)&(ebossALL['ra']<=28))]
	patches[4]=ebossALL[(ebossALL['ra']>28)]
	patches[5]=data[(data['programname']==b'boss')]
	patches[6]=data[(data['programname']==b'legacy')]

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
			dmm.append(None)
			rmm.append(None)

    #get areas
	area=np.empty(7, dtype=float)
	for i in range(7):
		if (len(patches[i]['z'])>1):
			if ((i==2)&(len(patches[1]['z'])>1)):
				area[2]=((rmm[2][1]-rmm[2][0])*(dmm[2][1]-dmm[2][0]))-((rmm[1][1]-rmm[2][0])*(dmm[1][1]-dmm[1][0]))
			else:
				area[i]=((rmm[i][1]-rmm[i][0])*(dmm[i][1]-dmm[i][0]))

	#calculate number of randoms to make for each patch
	ddensity=np.empty(8, dtype=float)
	temp=np.empty(8, dtype=float)
	pop=np.empty(8, dtype=int)
	for i, p in enumerate(patches): 
		if(len(p['z'])>1):
			ddensity[i]=dr*len(p)/area[i]
			pop[i]=dr*len(p)
			if i==2:
				rp2=pop[2]
				pop[2]=int(math.ceil(ddensity[2]*((rmm[2][1]-rmm[2][0])*(dmm[2][1]-dmm[2][0]))))

    #Create randoms for each patch
	random=[]
	for i, p in enumerate(patches):
		if(len(p['z'])>1):
		    r=[]
		    d=[]
		    z=[]
		    q=np.cumsum(evalu[i])/np.cumsum(evalu[i])[-1]
		    inv=interpolate.interp1d(q, b)
		    for j in range (pop[i]):
		        r.append(rand.uniform(rmm[i][0],rmm[i][1]))
		        d.append(rand.uniform(dmm[i][0],dmm[i][1]))
		        a=rand.uniform(q.min(),q.max())
		        z.append(inv(a))
		    temp=Table([r, d, z], names=('ra', 'dec', 'z'), meta={'name': 'first table'})
		    if i==2:
		        temp=temp[(((temp['ra']>0)&(temp['ra']<=14))|((((temp['ra']>=-3)&(temp['ra']<=0)))&((temp['dec']>1)|(temp['dec']<-1))))]
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

	temp=data[data['ra']<0]
	for t in temp:
		t['ra']+=360
	data[data['ra']<0]=temp

	total=np.array(total)
	total=z_to_cdist(total,cosmo=cosmo)

	return total




