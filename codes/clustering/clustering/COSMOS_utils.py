import numpy as np 
from numpy.lib.recfunctions import append_fields
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import SkyCoord
from astropy import units as u
import sys
from scipy import interpolate

path = '/Users/justin/Astro/Yale18/Clustering/clustering/COSMOS/'
fluxarea=np.genfromtxt(path+'flux_area_chandra_cosmos_legacy_05-10keV.dat')
F = fluxarea[:,0]
A = fluxarea[:,1]
f = interpolate.interp1d(F,A)
CosmosArea=2.161

def area(log_flux):
    if log_flux>-11:
        return 2.1609
    elif log_flux<-15.06:
        return 1e-6
    else:
        return f(log_flux)

def include_area_weights(cat):
    p_arr=[]
    for f in cat['flux']:
        p = area(np.log10(f))/CosmosArea
        p_arr.append(1./p)
    norm = len(cat)/np.sum(p_arr)
    w_arr=norm*np.array(p_arr)
    if 'weight' in cat.dtype.names:
        w_arr = w_arr * cat['weight']
    temp = list(zip(cat['z'],cat['ra'],cat['dec'],cat['flux'],w_arr))
    new = np.zeros((len(cat),), dtype=[('z', '<f8'),('ra', '<f8'),('dec', '<f8'),('flux', '<f8'),('weight', '<f8')])
    new[:] = temp
    return new


def parseCOSMOS_pdfs(data, pdf_dir=path+'chandra_cosmos_legacy_pdz/', extension='spec.light'):
    '''
    Input: data array with following column names:
    ra: 'CP_RA'
    dec: 'CP_DEC'
    spec-z: 'SPEC_Z'
    phot-z: 'PHOTO_Z'
    'PDZ'

    returns parsed array with z weights, array of PDFs
    '''
    if pdf_dir is None:
        sys.exit('Please specify directory for photo-z pdf files')

    w_arr=[]
    z_arr=[]
    ra_arr=[]
    dec_arr=[]
    id_arr=[]
    flux_arr=[]
    ztype=[]
    pdfs=[]
    skip=0
    for i,r in enumerate(data):
        idx = (r['id_x']).strip().decode()
        if (r['z_spec']<=0) and (r['z_phot']>0):
            fname = pdf_dir+idx+'.'+extension
            try:
                full_pdf = np.genfromtxt(fname)
            except IOError:
                print('unable to find PDF file')
                skip = skip+1
                continue

            pdf=full_pdf[np.where(full_pdf[:,1]>1e-4)]
            p=pdf[:,1]
            z=pdf[:,0]
            #if len(p)>100:
            #   p=p[::10]
            #   z=z[::10]
            p=p/sum(p)
            ra = r['ra_i']*np.ones(len(p))
            dec = r['dec_i']*np.ones(len(p))
            ident = np.chararray(len(p),32)
            ident[:] = idx
            #ident = idx*np.ones(len(p))
            flux = r['flux_F']*np.ones(len(p))
            phot = np.ones(len(p))
            w_arr=np.append(w_arr,p)
            z_arr=np.append(z_arr,z)
            ra_arr=np.append(ra_arr,ra)
            dec_arr=np.append(dec_arr,dec)
            flux_arr=np.append(flux_arr,flux)
            id_arr=np.append(id_arr,ident)
            ztype=np.append(ztype,phot)
            pdfs.append([z,p])
        else:
            w_arr=np.append(w_arr,1.)
            z_arr=np.append(z_arr,r['z_spec'])
            ra_arr=np.append(ra_arr,r['ra_i'])
            dec_arr=np.append(dec_arr,r['dec_i'])
            id_arr=np.append(id_arr,idx)
            flux_arr=np.append(flux_arr,r['flux_F'])
            ztype=np.append(ztype,0)
            pdfs.append([np.array([r['z_spec']]),np.array([1.])])
    temp = list(zip(ra_arr,dec_arr,z_arr,flux_arr,w_arr,id_arr,ztype))
    parsed = np.zeros((len(ra_arr),), dtype=[('ra', '<f8'),('dec', '<f8'),('z', '<f8'),('flux', '<f8'),('weight', '<f8'),('id', '<U32'),('ztype', '<i8')])
    parsed[:] = temp

    return parsed


