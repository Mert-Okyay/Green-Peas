import numpy as np
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table
from astropy.io import fits
import time
import sys
import corner
import emcee

from halotools.empirical_models import PrebuiltHodModelFactory, PrebuiltSubhaloModelFactory
from halotools.mock_observables import return_xyz_formatted_array
from halotools.mock_observables import wp
from halotools.sim_manager import CachedHaloCatalog

direc = '/home/mcp74/Science/Projects/clustering/bass/'

AG = np.array(Table.read(direc+'results/wp_ag.txt', format='ascii'))
rp = AG['rp']
wp_ag = AG['wp']
wp_ag_err = AG['wp_err']
wp_ag_cov = np.loadtxt(direc+'results/cov.txt')

rpbins = np.logspace(-1, np.log10(40), len(rp) + 1)

halocat = CachedHaloCatalog(simname = 'bolplanck', redshift = 0, halo_finder = 'rockstar')
nwalkers = 20
nsteps = 100

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

gall_positions_best = np.loadtxt(direc+'catalogs/gall_positions.txt')

halocat = CachedHaloCatalog(simname = 'bolplanck', redshift = 0, halo_finder = 'rockstar')

Amodel = PrebuiltHodModelFactory('zheng07',redshift=0,threshold=-18)
const = Amodel.param_dict['logM1'] -Amodel.param_dict['logM0']
Amodel.populate_mock(halocat)

def wp_HOD_x(rpbins,Mcmin,alpha):
    pimax=60
    
    Amodel.param_dict['sigma_logM']=0
    Amodel.param_dict['alpha']=alpha
    Amodel.param_dict['logM0']=Mcmin
    Amodel.param_dict['logM1']=Mcmin+const
    Amodel.param_dict['logMmin']=Mcmin
    #print(alpha,Mmin,M0)

    #Amodel.populate_mock(halocat)
    Amodel.mock.populate()
    
    x = Amodel.mock.galaxy_table['x']
    y = Amodel.mock.galaxy_table['y']
    z = Amodel.mock.galaxy_table['z']
    all_positions = return_xyz_formatted_array(x, y, z)

    #print len(x)
    w = wp(all_positions, rpbins, pimax, sample2=gall_positions_best, period=Amodel.mock.Lbox, \
           num_threads=4, do_cross=True)
    return w[1]

def lnprior(theta):
    mmin,alpha = theta
    if 11 < mmin < 14 and -.5 < alpha < 1.5:
        return 0.0
    return -np.inf

invcov = np.linalg.inv(wp_ag_cov)

def lnprob(theta, rpbins, wpd, icov):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    else:
        Mcmin,alpha=theta
        wpm = wp_HOD_x(rpbins,Mcmin,alpha)
        diff = wpm-wpd
        return -np.dot(diff,np.dot(icov,diff))/2.0


ndim = 2
p0 = np.loadtxt(direc+'mcmc/pos.txt')
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(rpbins, wp_ag, invcov))

f = open(direc+"chain.dat", "w")
f.close()

for result in sampler.sample(p0, iterations=500, storechain=False):
    position = result[0]
    f = open(direc+"chain.dat", "a")
    for k in range(position.shape[0]):
        f.write("{0:4d} {1:s}\n".format(k, " ".join(str(position[k]))))
    f.close()