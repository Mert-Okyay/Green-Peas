
from __future__ import division

import numpy as np 
from numpy.lib.recfunctions import append_fields
from astropy.table import Table, vstack
from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import SkyCoord
from astropy import units as u
import sys
from joblib import Parallel, delayed
import math
import matplotlib.pyplot as plt
from clustering.clustering.utils import *
from Corrfunc.mocks.DDrppi_mocks import DDrppi_mocks

def auto_jackknife(d,r,m,pimax,bins,estimator,covariance=True,survey=None, JKN=True):
	#print('beginning jackknifes')
	if survey=='BASS':
		dblocks=_BASS_block(d,m)
		rblocks=_BASS_block(r,m)
	elif survey=='S82X':
		dblocks=_S82X_block(d,m)
		rblocks=_S82X_block(r,m)
	elif survey=='S82X_test1':
		dblocks=_S82X_block_test1(d,m)
		rblocks=_S82X_block_test1(r,m)
	elif survey=='S82X_test2':
		dblocks=_S82X_block_test2(d,m)
		rblocks=_S82X_block_test2(r,m)
	elif survey=='S82X_test3':
		dblocks=_S82X_block_test3(d,m)
		rblocks=_S82X_block_test3(r,m)
	elif survey=='AO13':
		dblocks=_AO13_block(d,m)
		rblocks=_AO13_block(r,m)
	elif survey=='COSMOS':
		dblocks=_COSMOS_block(d,m)
		rblocks=_COSMOS_block(r,m)
	elif survey=='XMM_XXL':
		dblocks=_XMM_XXL_block(d,m)
		rblocks=_XMM_XXL_block(r,m)
	elif survey=='Comb':
		dblocks=_Comb_block(d,m)
		rblocks=_Comb_block(r,m)
	elif survey=='zCosDEEP':
		dblocks=_zCosDEEP_block(d,m)
		rblocks=_zCosDEEP_block(r,m)
	else:
		print('No valid survey for error estimation')
		return
	M=len(dblocks)
	if JKN:
		print('Using ',M,' jacknife samples')
	nbins = len(bins)-1
	
	#NON-PARALLELIZED JACKKNIFE:
	wp_arr = np.empty(nbins)
	for i in np.arange(M):
		if JKN:
			print(str(i+1),'/',str(M))
		new_dblocks = np.delete(dblocks,i,0)
		new_rblocks = np.delete(rblocks,i,0)
		dset = np.concatenate(new_dblocks)
		rset = np.concatenate(new_rblocks)
		wp_t = wp_dd(data=dset,randoms=rset,pimax=pimax,bins=bins,estimator=estimator)
		wp_arr=np.vstack((wp_arr,wp_t))
	wp_arr=wp_arr[1:]
	
	#PARALLELIZED JACKKNIFE COMPUTATION:
	#wp_arr = np.array(Parallel(n_jobs=-1)(delayed(PAjackknife)(i,dblocks,rblocks,pimax=pimax,bins=bins,estimator=estimator) for i in np.arange(M)))

	err =np.sqrt(((M-1.)/M) * np.sum( (wp_arr - np.average(wp_arr,axis=0))**2, axis=0))
	if covariance:
		cov = np.zeros(shape = (nbins,nbins))
		for i in np.arange(nbins):
			for j in np.arange(nbins):
				cov[i,j] = ((M-1.)/M)* np.sum( (wp_arr[:,i] - np.average(wp_arr,axis=0)[i]) * (wp_arr[:,j] - np.average(wp_arr,axis=0)[j]), axis=0)
		return err,cov
	else: return err

def cross_jackknife(d1,d2,r1,r2,m,pimax,bins,estimator,covariance=True,survey=None,weights1=None,weights2=None, adaptive=0, JKN=True):
	#print('beginning jackknifes')
	if survey=='BASS':
		d1blocks=_BASS_block(d1,m)
		d2blocks=_BASS_block(d2,m)
		r1blocks=_BASS_block(r1,m)
		r2blocks=_BASS_block(r2,m)
	elif survey=='S82X':
		d1blocks=_S82X_block(d1,m)
		d2blocks=_S82X_block(d2,m)
		r1blocks=_S82X_block(r1,m)
		r2blocks=_S82X_block(r2,m)
	elif survey=='S82X_test1':
		d1blocks=_S82X_block_test1(d1,m)
		d2blocks=_S82X_block_test1(d2,m)
		r1blocks=_S82X_block_test1(r1,m)
		r2blocks=_S82X_block_test1(r2,m)
	elif survey=='S82X_test2':
		d1blocks=_S82X_block_test2(d1,m)
		d2blocks=_S82X_block_test2(d2,m)
		r1blocks=_S82X_block_test2(r1,m)
		r2blocks=_S82X_block_test2(r2,m)
	elif survey=='S82X_test3':
		d1blocks=_S82X_block_test3(d1,m)
		d2blocks=_S82X_block_test3(d2,m)
		r1blocks=_S82X_block_test3(r1,m)
		r2blocks=_S82X_block_test3(r2,m)
	elif survey=='AO13':
		d1blocks=_AO13_block(d1,m)
		d2blocks=_AO13_block(d2,m)
		r1blocks=_AO13_block(r1,m)
		r2blocks=_AO13_block(r2,m)
	elif survey=='COSMOS':
		d1blocks=_COSMOS_block(d1,m)
		d2blocks=_COSMOS_block(d2,m)
		r1blocks=_COSMOS_block(r1,m)
		r2blocks=_COSMOS_block(r2,m)
	elif survey=='XMM_XXL':
		d1blocks=_XMM_XXL_block(d1,m)
		d2blocks=_XMM_XXL_block(d2,m)
		r1blocks=_XMM_XXL_block(r1,m)
		r2blocks=_XMM_XXL_block(r2,m)
	elif survey=='Comb':
		d1blocks=_Comb_block(d1,m)
		d2blocks=_Comb_block(d2,m)
		r1blocks=_Comb_block(r1,m)
		r2blocks=_Comb_block(r2,m)
	elif survey=='zCosDEEP':
		d1blocks=_zCosDEEP_block(d1,m)
		d2blocks=_zCosDEEP_block(d2,m)
		r1blocks=_zCosDEEP_block(r1,m)
		r2blocks=_zCosDEEP_block(r2,m)
	else:
		print('No valid survey for error estimation')
		return
	M=len(d1blocks)
	if(JKN):
		print('Using ',M,' jacknife samples')
	nbins = len(bins)-1
	
	
	#compute wp without one block
	wp_arr = np.empty(nbins)
	for i in np.arange(M):
		if(JKN):
			print(str(i+1),'/',str(M))
		new_d1blocks = np.delete(d1blocks,i,0)
		new_d2blocks = np.delete(d2blocks,i,0)
		if r1 is not None:
			new_r1blocks = np.delete(r1blocks,i,0)
		new_r2blocks = np.delete(r2blocks,i,0)
		d1set = np.concatenate(new_d1blocks)
		d2set = np.concatenate(new_d2blocks)
		if r1 is not None:
			r1set = np.concatenate(new_r1blocks)
		else:
			r1set=None
		r2set = np.concatenate(new_r2blocks)
		wp_t, binst = wp_d1d2(d1=d1set, d2=d2set, r1=r1set, r2=r2set, bins=bins, pimax=pimax, estimator=estimator)
		wp_arr=np.vstack((wp_arr,wp_t))
	wp_arr=wp_arr[1:]

	# For parallel processing of Jackknife samples:
	#wp_arr = np.array(Parallel(n_jobs=-1)(delayed(PCjackknife)(i,d1blocks,d2blocks,r2blocks,pimax=pimax,bins=bins,estimator=estimator,r1blocks=r1blocks) for i in np.arange(M)))

	# For non-parallel processing of Jackknife sample:
	err =np.sqrt(((M-1.)/M) * np.sum( (wp_arr - np.average(wp_arr,axis=0))**2, axis=0))
	if covariance:
		cov = np.zeros(shape = (nbins,nbins))
		for i in np.arange(nbins):
			for j in np.arange(nbins):
				cov[i,j] = ((M-1.)/M)* np.sum( (wp_arr[:,i] - np.average(wp_arr,axis=0)[i]) * (wp_arr[:,j] - np.average(wp_arr,axis=0)[j]), axis=0)
		return err,cov
	else: return err

def PAjackknife(i,dblocks,rblocks,pimax,bins,estimator):
	new_dblocks = np.delete(dblocks,i,0)
	new_rblocks = np.delete(rblocks,i,0)
	dset = np.concatenate(new_dblocks)
	rset = np.concatenate(new_rblocks)
	wp_t = wp_dd(data=dset,randoms=rset,pimax=pimax,bins=bins,estimator=estimator)
	return wp_t

def PCjackknife(i,d1blocks,d2blocks,r2blocks,pimax,bins,estimator,r1blocks=None):
	new_d1blocks = np.delete(d1blocks,i,0)
	new_d2blocks = np.delete(d2blocks,i,0)
	#print(r1blocks)
	if r1blocks is not 0:
		new_r1blocks = np.delete(r1blocks,i,0)
	new_r2blocks = np.delete(r2blocks,i,0)
	d1set = np.concatenate(new_d1blocks)
	d2set = np.concatenate(new_d2blocks)
	if r1blocks is not 0:
		r1set = np.concatenate(new_r1blocks)
	else:
		r1set=None
	r2set = np.concatenate(new_r2blocks)
	wp_t = wp_d1d2(d1=d1set, d2=d2set, r1=r1set, r2=r2set, bins=bins, pimax=pimax, estimator=estimator)
	return wp_t

def _BASS_block(cat,m):

	if cat==None:
		return 0
	blocks=[]
	lens=[]

	ra0=0
	ra_int=(360./m)
	dec_int=(180./m)
	#dec_int=(2./m)
	for i in range(m):
		dec0=-90
		#dec0=-1
		for j in range(m):
			ra_min=ra0
			ra_max=ra0+ra_int
			dec_min=dec0
			dec_max=dec0+dec_int
			#dec_min = (np.arcsin(dec0)*u.radian).to(u.degree).value
			#dec_max = (np.arcsin(dec0+dec_int)*u.radian).to(u.degree).value
			block=cat[(cat['ra']<ra_max)&(cat['ra']>=ra_min)&\
					(cat['dec']<dec_max)&(cat['dec']>=dec_min)]
			dec0=dec_max
			blocks.append(block)
			lens.append(len(block))
			#print ra_min,dec_min
		ra0=ra_max
	return blocks

def _AO13_block(cat,m):
	#if cat==None:
	#	return 0
	#cat.sort(order='ra')
	#x = len(cat)%m
	#if x!=0: newc = cat[:-x]
	#else: newc=cat
	#blocks = np.array(np.split(newc, m))
	#return blocks

	if cat is None:
		return 0
	blocks=[]
	lens=[]

	ramin=14.1
	ramax=28.04
	ra_int=((ramax-ramin)/m)
	ra0=ramin
	dec_min=-0.63
	dec_max=0.63
	for i in range(m):
		ra_min=ra0
		ra_max=ra0+ra_int
		block=cat[(cat['ra']<ra_max)&(cat['ra']>=ra_min)&\
				(cat['dec']<dec_max)&(cat['dec']>=dec_min)]
		dec0=dec_max
		blocks.append(block)
		lens.append(len(block))
		#print ra_min,dec_min
		ra0=ra_max
	return blocks

def _zCosDEEP_block(cat,m):
	if cat is None:
		return 0
	blocks=[]
	lens=[]

	ramax = 150.645706
	ramin = 149.710846
	dec_min = 1.718882
	dec_max = 2.63344

	ra_int=((ramax-ramin)/m)
	ra0=ramin

	for i in range(m):
		ra_min=ra0
		ra_max=ra0+ra_int
		block=cat[(cat['ra']<ra_max)&(cat['ra']>=ra_min)&\
				(cat['dec']<dec_max)&(cat['dec']>=dec_min)]
		dec0=dec_max
		blocks.append(block)
		lens.append(len(block))
		ra0=ra_max

	# print('percentages:')
	# temp=[]
	# for i, b in enumerate(blocks):
	# 	print('block %2d: %f' %(i, len(b)/len(cat)*100))
	# 	temp.append[len(b)/len(cat)*100]
	# print('range: %f\t%f' %(temp.min(),temp.max()))
	return blocks

def _S82X_block(cat, m):
	if cat is None:
		return 0
	blocks=[]
	blocks.append(cat[(cat['ra']>300)&(cat['ra']<=332)])
	blocks.append(cat[(cat['ra']>332)&(cat['ra']<=334.3152435687774)])
	blocks.append(cat[(cat['ra']>334.3152435687774)&(cat['ra']<=340)])
	blocks.append(cat[(cat['ra']>340)&(cat['ra']<=352.50688215)])
	blocks.append(cat[(cat['ra']>352.50688215)&(cat['ra']<=353.5)])
	blocks.append(cat[(cat['ra']>353.5)&(cat['ra']<=358)])
	blocks.append(cat[(cat['ra']>358)|(cat['ra']<=4)])
	blocks.append(cat[(cat['ra']>4)&(cat['ra']<14.1)])
	ramin=14.1
	ramax=28.04
	ra_int=((ramax-ramin)/10)
	ra0=ramin
	dec_min=-0.63
	dec_max=0.63
	for i in range(10):
		ra_min=ra0
		ra_max=ra0+ra_int
		block=cat[(cat['ra']<ra_max)&(cat['ra']>=ra_min)]
		dec0=dec_max
		blocks.append(block)
		ra0=ra_max
	blocks.append(cat[(cat['ra']>28.04)&(cat['ra']<=42)])
	blocks.append(cat[(cat['ra']>42)&(cat['ra']<=50)])
	blocks.append(cat[(cat['ra']>50)&(cat['ra']<=60)])
	return blocks

def _S82X_block_test1(cat, m):
	if cat is None:
		return 0
	blocks=[]
	blocks.append(cat[(cat['ra']>300)&(cat['ra']<=334.3152435687774)])
	blocks.append(cat[(cat['ra']>334.3152435687774)&(cat['ra']<=352.50688215)])
	blocks.append(cat[(cat['ra']>352.50688215)&(cat['ra']<=358)])
	blocks.append(cat[(cat['ra']>358)|(cat['ra']<14.1)])
	ramin=14.1
	ramax=28.04
	ra_int=((ramax-ramin)/5)
	ra0=ramin
	dec_min=-0.63
	dec_max=0.63
	for i in range(5):
		ra_min=ra0
		ra_max=ra0+ra_int
		block=cat[(cat['ra']<ra_max)&(cat['ra']>=ra_min)]
		dec0=dec_max
		blocks.append(block)
		ra0=ra_max
	blocks.append(cat[(cat['ra']>28.04)&(cat['ra']<=50)])
	blocks.append(cat[(cat['ra']>50)&(cat['ra']<=60)])
	return blocks

def _S82X_block_test2(cat, m):
	if cat is None:
		return 0
	blocks=[]
	blocks.append(cat[(cat['ra']>300)&(cat['ra']<=332)])
	blocks.append(cat[(cat['ra']>332)&(cat['ra']<=340)])
	blocks.append(cat[(cat['ra']>340)&(cat['ra']<=353.5)])
	blocks.append(cat[(cat['ra']>353.5)|(cat['ra']<=4)])
	blocks.append(cat[(cat['ra']>4)&(cat['ra']<14.1)])
	ramin=14.1
	ramax=28.04
	ra_int=((ramax-ramin)/5)
	ra0=ramin
	for i in range(5):
		ra_min=ra0
		ra_max=ra0+ra_int
		block=cat[(cat['ra']<ra_max)&(cat['ra']>=ra_min)]
		blocks.append(block)
		ra0=ra_max
	blocks.append(cat[(cat['ra']>28.04)&(cat['ra']<=42)])
	blocks.append(cat[(cat['ra']>42)&(cat['ra']<=60)])
	return blocks

def _S82X_block_test3(dsp2, m):
	if dsp2 is None:
		return 0
	patch=[]
	patch.append(dsp2[(dsp2['ra']>300)&(dsp2['ra']<311)])
	patch.append(dsp2[(dsp2['ra']>311)&(dsp2['ra']<323)])
	patch.append(dsp2[(dsp2['ra']>323)&(dsp2['ra']<324)])
	patch.append(dsp2[(dsp2['ra']>324)&(dsp2['ra']<332)])
	patch.append(dsp2[(dsp2['ra']>332)&(dsp2['ra']<334.3152435687774)])
	patch.append(dsp2[(dsp2['ra']>334.3152435687774)&(dsp2['ra']<336)])#bps
	patch.append(dsp2[(dsp2['ra']>336)&(dsp2['ra']<350)])
	patch.append(dsp2[(dsp2['ra']>350)&(dsp2['ra']<354)])#bps
	patch.append(dsp2[(dsp2['ra']>354)&(dsp2['ra']<355)])
	patch.append(dsp2[(dsp2['ra']>355)&(dsp2['ra']<356.5)])
	patch.append(dsp2[(dsp2['ra']>356.5)&(dsp2['ra']<358)])
	patch.append(dsp2[(dsp2['ra']>358)&(dsp2['ra']<360)])
	patch.append(dsp2[(dsp2['ra']>0)&(dsp2['ra']<2)])
	patch.append(dsp2[(dsp2['ra']>2)&(dsp2['ra']<4)])
	patch.append(dsp2[(dsp2['ra']>4)&(dsp2['ra']<7)])
	patch.append(dsp2[(dsp2['ra']>7)&(dsp2['ra']<10.3)])
	patch.append(dsp2[(dsp2['ra']>10.3)&(dsp2['ra']<14.1)])
	ramin=14.1
	ramax=28.04
	ra_int=((ramax-ramin)/20)
	ra0=ramin
	dec_min=-0.63
	dec_max=0.63
	for i in range(20):
		ra_min=ra0
		ra_max=ra0+ra_int
		block=dsp2[(dsp2['ra']<ra_max)&(dsp2['ra']>=ra_min)]
		dec0=dec_max
		patch.append(block)
		ra0=ra_max
	patch.append(dsp2[(dsp2['ra']>28.04)&(dsp2['ra']<29)])
	patch.append(dsp2[(dsp2['ra']>29)&(dsp2['ra']<34)])
	patch.append(dsp2[(dsp2['ra']>34)&(dsp2['ra']<39)])
	patch.append(dsp2[(dsp2['ra']>39)&(dsp2['ra']<42)])
	patch.append(dsp2[(dsp2['ra']>42)&(dsp2['ra']<45)])
	patch.append(dsp2[(dsp2['ra']>45)&(dsp2['ra']<46)&(dsp2['dec']>-0.6)])
	patch.append(dsp2[(dsp2['ra']>46)&(dsp2['ra']<49)&(dsp2['dec']>-0.6)])
	patch.append(dsp2[(dsp2['ra']>45)&(dsp2['ra']<49)&(dsp2['dec']<-0.6)])
	patch.append(dsp2[(dsp2['ra']>49)&(dsp2['ra']<53.9)])
	patch.append(dsp2[(dsp2['ra']>53.9)&(dsp2['ra']<54.38)])
	patch.append(dsp2[(dsp2['ra']>54.38)&(dsp2['ra']<55)])
	patch.append(dsp2[(dsp2['ra']>55)&(dsp2['ra']<56)])
	patch.append(dsp2[(dsp2['ra']>56)&(dsp2['ra']<60)])
	return patch

def _COSMOS_block(cat, m):
	if cat is None:
		return 0
	print_perct=False
	blocks=[]
	ramin=149.319455
	ramax=150.918003
	decmin=1.419393
	decmax=3.0089
	n=math.ceil(math.sqrt(m))
	ra_int=((ramax-ramin)/n)
	dec_int=((decmax-decmin)/n)
	ra0=ramin
	dec0=decmin
	for i in range(n):
		dec_min=dec0
		dec_max=dec0+dec_int
		for j in range(n):
			ra_min=ra0
			ra_max=ra0+ra_int
			block=cat[(cat['ra']<ra_max)&(cat['ra']>=ra_min)&(cat['dec']<dec_max)&(cat['dec']>=dec_min)]
			blocks.append(block)
			ra0=ra_max
		ra0=ramin
		dec0=dec_max
	if print_perct:
		for i, b in enumerate(blocks):
			if 'weight' in cat.dtype.names:
				print(i, (sum(b['weight'])/sum(cat['weight']))*100)
			else:
				print(i, (len(b)/len(cat))*100)
	return blocks

def _XMM_XXL_grid_block(cat, m):
	printcat=False
	if cat is None:
		return 0
	blocks=[]
	ramin=30.0592
	ramax=38.8278
	decmin=-7.34838
	decmax=-2.67178
	n=math.ceil(math.sqrt(m))
	ra_int=((ramax-ramin)/n)
	dec_int=((decmax-decmin)/n)
	dec0=decmin

	if printcat:
		fig, ax = plt.subplots()
		gra = np.array(cat['ra'])
		gdec = np.array(cat['dec'])
		ax.scatter(gra,gdec,1.,color='k',alpha=.5)
		#ax.axis([30.0593, 38.8277, -7.34839, -2.67177])
		ax.set_xlabel('ra')
		ax.set_ylabel('dec')

	for i in range(n):
		dec_min=dec0
		dec_max=dec0+dec_int
		ra0=ramin
		for j in range(n):
			ra_min=ra0
			ra_max=ra0+ra_int
			block=cat[(cat['ra']<ra_max)&(cat['ra']>=ra_min)&(cat['dec']<dec_max)&(cat['dec']>=dec_min)]
			if printcat:
				gra = np.array(block['ra'])
				gdec = np.array(block['dec'])
				ax.scatter(gra,gdec,1.,alpha=.5)
			blocks.append(block)
			ra0=ra_max
		dec0=dec_max

	return blocks

def _XMM_XXL_block(cat, m):
	if cat is None:
		return 0
	blocks=[]
	ramin=30.0592
	ramax=38.8278
	decmin=-7.34838
	decmax=-2.67178
	#n=math.ceil(math.sqrt(m))
	n=5
	ra_int=((ramax-ramin)/n)
	dec_int=((decmax-decmin)/n)
	dec0=decmin
	patch=[]
	for i in range(n):
		dec_min=dec0
		dec_max=dec0+dec_int
		ra0=ramin
		for j in range(n):
			ra_min=ra0
			ra_max=ra0+ra_int
			new=cat[(cat['ra']<ra_max)&(cat['ra']>=ra_min)&(cat['dec']<dec_max)&(cat['dec']>=dec_min)]
			if len(patch)==0:
				patch=new
			else:
				patch=np.array(vstack([Table(patch),Table(new)]))
			if (((i==1)&(j==1))|((i==1)&(j==4))|((i==2)&(j>0))|((i==3)&(j>0)&(j<4))|((i==4)&(j==3))):
				blocks.append(patch)
				patch=[]
			elif ((i==3)&(j==4)):
				side=patch
				patch=[]
			elif((i==4)&(j==4)):
				patch=np.array(vstack([Table(patch),Table(side)]))
				blocks.append(patch)
				patch=[]
			ra0=ra_max
		dec0=dec_max
	blocks=np.array(blocks)
	return blocks

def _Comb_block(cat, m):
	if cat is None:
		return 0
	perct=False
	blocks=[]
	#Stripe82

	blocks.append(cat[(cat['ra']>300)&(cat['ra']<=332)&(cat['dec']<=1.5)&(cat['dec']>=-1.5)])
	blocks.append(cat[(cat['ra']>332)&(cat['ra']<=334.3152435687774)&(cat['dec']<=1.5)&(cat['dec']>=-1.5)])
	blocks.append(cat[(cat['ra']>334.3152435687774)&(cat['ra']<=340)&(cat['dec']<=1.5)&(cat['dec']>=-1.5)])
	blocks.append(cat[(cat['ra']>340)&(cat['ra']<=352.50688215)&(cat['dec']<=1.5)&(cat['dec']>=-1.5)])
	blocks.append(cat[(cat['ra']>352.50688215)&(cat['ra']<=353.5)&(cat['dec']<=1.5)&(cat['dec']>=-1.5)])
	blocks.append(cat[(cat['ra']>353.5)&(cat['ra']<=358)&(cat['dec']<=1.5)&(cat['dec']>=-1.5)])
	blocks.append(cat[((cat['ra']>358)|(cat['ra']<=4))&(cat['dec']<=1.5)&(cat['dec']>=-1.5)])
	blocks.append(cat[(cat['ra']>4)&(cat['ra']<14.1)&(cat['dec']<=1.5)&(cat['dec']>=-1.5)])
	ramin=14.1
	ramax=28.04
	ra_int=((ramax-ramin)/10)
	ra0=ramin
	for i in range(10):
		ra_min=ra0
		ra_max=ra0+ra_int
		block=cat[(cat['ra']<ra_max)&(cat['ra']>=ra_min)&(cat['dec']<=1.5)&(cat['dec']>=-1.5)]
		blocks.append(block)
		ra0=ra_max
	blocks.append(cat[(cat['ra']>28.04)&(cat['ra']<=42)&(cat['dec']<=1.5)&(cat['dec']>=-1.5)])
	blocks.append(cat[(cat['ra']>42)&(cat['ra']<=50)&(cat['dec']<=1.5)&(cat['dec']>=-1.5)])
	blocks.append(cat[(cat['ra']>50)&(cat['ra']<=60)&(cat['dec']<=1.5)&(cat['dec']>=-1.5)])
	#XMM-XXL
	ramin=30.0592
	ramax=38.8278
	decmin=-7.34838
	decmax=-2.67178
	#n=math.ceil(math.sqrt(m))
	n=5
	ra_int=((ramax-ramin)/n)
	dec_int=((decmax-decmin)/n)
	dec0=decmin
	patch=[]
	for i in range(n):
		dec_min=dec0
		dec_max=dec0+dec_int
		ra0=ramin
		for j in range(n):
			ra_min=ra0
			ra_max=ra0+ra_int
			new=cat[(cat['ra']<ra_max)&(cat['ra']>=ra_min)&(cat['dec']<dec_max)&(cat['dec']>=dec_min)]
			if len(patch)==0:
				patch=new
			else:
				patch=np.array(vstack([Table(patch),Table(new)]))
			if (((i==1)&(j==1))|((i==1)&(j==4))|((i==2)&(j>0))|((i==3)&(j>0)&(j<4))|((i==4)&(j==3))):
				blocks.append(patch)
				patch=[]
			elif ((i==3)&(j==4)):
				side=patch
				patch=[]
			elif((i==4)&(j==4)):
				patch=np.array(vstack([Table(patch),Table(side)]))
				blocks.append(patch)
				patch=[]
			ra0=ra_max
		dec0=dec_max
	#COSMOS

	# ramin=149.319455
	# ramax=150.918003
	# decmin=1.419393
	# decmax=3.0089
	# n=4
	# ra_int=((ramax-ramin)/n)
	# dec_int=((decmax-decmin)/n)
	# ra0=ramin
	# dec0=decmin
	# for i in range(n):
	# 	dec_min=dec0
	# 	dec_max=dec0+dec_int
	# 	for j in range(n):
	# 		ra_min=ra0
	# 		ra_max=ra0+ra_int
	# 		block=cat[(cat['ra']<ra_max)&(cat['ra']>=ra_min)&(cat['dec']<dec_max)&(cat['dec']>=dec_min)]
	# 		blocks.append(block)
	# 		ra0=ra_max
	# 	ra0=ramin 
	# 	dec0=dec_max
	if perct:
		for i, b in enumerate(blocks):
			if 'weight' in cat.dtype.names:
				print(i, (sum(b['weight'])/sum(cat['weight']))*100)
			else:
				print(i, (len(b)/len(cat))*100)
	return blocks





