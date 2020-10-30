import numpy as np 
from numpy.lib.recfunctions import append_fields
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import SkyCoord
from astropy import units as u
import sys
from scipy import interpolate

path = '/Users/justin/Astro/Yale18/Clustering/clustering/stripe82x/'
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


def parseS82X_pdfs(data, pdf_dir, prefix='Id', extension='spec', zfill=9):
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
    if 'nh' in data.dtype.names:
        w_arr=[]
        z_arr=[]
        ra_arr=[]
        dec_arr=[]
        id_arr=[]
        flux_arr=[]
        hr_arr=[]
        ztype=[]
        pdfs=[]
        nh_arr=[]
        skip=0
        for i,r in enumerate(data):
            if not (r['HARD_COUNTS']+r['SOFT_COUNTS'])==0:
                if (r['SPEC_Z']<=0) and (r['PHOTO_Z']>0):
                    fname = pdf_dir+prefix+str(r['REC_NO']).zfill(zfill)+'.'+extension
                    try:
                        full_pdf = np.genfromtxt(fname,skip_header=32)[:651]
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
                    ra = r['CP_RA']*np.ones(len(p))
                    dec = r['CP_DEC']*np.ones(len(p))
                    ident = r['REC_NO']*np.ones(len(p))
                    flux = r['FULL_FLUX']*np.ones(len(p))
                    nh = r['nh']*np.ones(len(p))

                    hr=((r['HARD_COUNTS']-r['SOFT_COUNTS'])/(r['HARD_COUNTS']+r['SOFT_COUNTS']))*np.ones(len(p))
                    hr_arr=np.append(hr_arr, hr)

                    phot = np.ones(len(p))
                    w_arr=np.append(w_arr,p)
                    z_arr=np.append(z_arr,z)
                    ra_arr=np.append(ra_arr,ra)
                    dec_arr=np.append(dec_arr,dec)
                    flux_arr=np.append(flux_arr,flux)
                    hr_arr=np.append(hr_arr, hr)
                    id_arr=np.append(id_arr,ident)
                    ztype=np.append(ztype,phot)
                    nh_arr=np.append(nh_arr,nh)
                    #Because dz is 0.02 for z>6:
                    #lz = (z<=6)
                    #hz = (z>6)
                    #new1=p[lz]
                    #new2=.5*np.repeat(p[hz],2.)
                    #newp=np.append(new1,new2)
                    #pdfs.append(newp)
                    pdfs.append([z,p])
                else:
                    w_arr=np.append(w_arr,1.)
                    z_arr=np.append(z_arr,r['SPEC_Z'])
                    ra_arr=np.append(ra_arr,r['CP_RA'])
                    dec_arr=np.append(dec_arr,r['CP_DEC'])
                    id_arr=np.append(id_arr,r['REC_NO'])
                    flux_arr=np.append(flux_arr,r['FULL_FLUX'])

                    hr=(r['HARD_COUNTS']-r['SOFT_COUNTS'])/(r['HARD_COUNTS']+r['SOFT_COUNTS'])
                    hr_arr=np.append(hr_arr, hr)

                    ztype=np.append(ztype,0)
                    nh_arr=np.append(nh_arr,r['nh'])
                    pdfs.append([np.array([r['SPEC_Z']]),np.array([1.])])
        if 'weight' in data.dtype.names:
            w_arr = w_arr * data['weight']
        temp = list(zip(ra_arr,dec_arr,z_arr,flux_arr,w_arr,id_arr,ztype, hr_arr,nh_arr))
        parsed = np.zeros((len(ra_arr),), dtype=[('ra', '<f8'),('dec', '<f8'),('z', '<f8'),('flux', '<f8'),('weight', '<f8'),('id', '<f8'),('ztype', '<i8'), ('hr', '<f8'), ('nh', '<f8')])
        parsed[:] = temp
    else:
        w_arr=[]
        z_arr=[]
        ra_arr=[]
        dec_arr=[]
        id_arr=[]
        flux_arr=[]
        hr_arr=[]
        ztype=[]
        pdfs=[]
        nh_arr=[]
        skip=0
        for i,r in enumerate(data):
            if not (r['HARD_COUNTS']+r['SOFT_COUNTS'])==0:
                if (r['SPEC_Z']<=0) and (r['PHOTO_Z']>0):
                    fname = pdf_dir+prefix+str(r['REC_NO']).zfill(zfill)+'.'+extension
                    try:
                        full_pdf = np.genfromtxt(fname,skip_header=32)[:651]
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
                    ra = r['CP_RA']*np.ones(len(p))
                    dec = r['CP_DEC']*np.ones(len(p))
                    ident = r['REC_NO']*np.ones(len(p))
                    flux = r['FULL_FLUX']*np.ones(len(p))
                    nh = -99*np.ones(len(p))

                    hr=((r['HARD_COUNTS']-r['SOFT_COUNTS'])/(r['HARD_COUNTS']+r['SOFT_COUNTS']))*np.ones(len(p))
                    hr_arr=np.append(hr_arr, hr)

                    phot = np.ones(len(p))
                    w_arr=np.append(w_arr,p)
                    z_arr=np.append(z_arr,z)
                    ra_arr=np.append(ra_arr,ra)
                    dec_arr=np.append(dec_arr,dec)
                    flux_arr=np.append(flux_arr,flux)
                    hr_arr=np.append(hr_arr, hr)
                    id_arr=np.append(id_arr,ident)
                    ztype=np.append(ztype,phot)
                    nh_arr=np.append(nh_arr,nh)
                    #Because dz is 0.02 for z>6:
                    #lz = (z<=6)
                    #hz = (z>6)
                    #new1=p[lz]
                    #new2=.5*np.repeat(p[hz],2.)
                    #newp=np.append(new1,new2)
                    #pdfs.append(newp)
                    pdfs.append([z,p])
                else:
                    w_arr=np.append(w_arr,1.)
                    z_arr=np.append(z_arr,r['SPEC_Z'])
                    ra_arr=np.append(ra_arr,r['CP_RA'])
                    dec_arr=np.append(dec_arr,r['CP_DEC'])
                    id_arr=np.append(id_arr,r['REC_NO'])
                    flux_arr=np.append(flux_arr,r['FULL_FLUX'])

                    hr=(r['HARD_COUNTS']-r['SOFT_COUNTS'])/(r['HARD_COUNTS']+r['SOFT_COUNTS'])
                    hr_arr=np.append(hr_arr, hr)

                    ztype=np.append(ztype,0)
                    nh_arr=np.append(nh_arr,-99)
                    pdfs.append([np.array([r['SPEC_Z']]),np.array([1.])])
        if 'weight' in data.dtype.names:
            w_arr = w_arr * data['weight']
        temp = list(zip(ra_arr,dec_arr,z_arr,flux_arr,w_arr,id_arr,ztype, hr_arr,nh_arr))
        parsed = np.zeros((len(ra_arr),), dtype=[('ra', '<f8'),('dec', '<f8'),('z', '<f8'),('flux', '<f8'),('weight', '<f8'),('id', '<f8'),('ztype', '<i8'), ('hr', '<f8'), ('nh', '<f8')])
        parsed[:] = temp

    return parsed, np.array(pdfs)

