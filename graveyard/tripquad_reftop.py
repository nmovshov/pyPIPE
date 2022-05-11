################################################################################
# Template module. Copy and modify before running.
# For MPI use pooltype='mpi' and run with:
# mpirun -n N python filename.py
# Remember: one mpi process will be master, N-1 will run lnprob
################################################################################
from __future__ import division
import os, sys
import pickle
import numpy as np
import emcee
import schwimmbad
import tof4
import dprofs
import losses
import priors
from timeit import default_timer as timer
from warnings import filterwarnings
#filterwarnings('ignore', category=RuntimeWarning, message='invalid value')
cout = sys.stdout.write
cerr = sys.stderr.write

### MODIFY LINES IN THIS SECTION ###
obsfile = r'Saturn_Obs.pickle'
seedfile = r'tqref3-mank1886-seed.txt'
reffile = r'tqref3-mank1886-ref.txt'
seedind = 0
mdl = dprofs.tripquad_reftop
lfh = losses.euclid_J24
pfh = priors.tripquad_reftop
toflevels = 2048
nsamples = 2
batchsize = 1
nwalkers = 20
ascale = 1.5
sskip = 32
pooltype = 'mpi'
fakelike = False
###

# Load reference top atmosphere
ref = np.loadtxt(reffile)
assert ref[-1,0]/ref[0,0] > 1/toflevels

# Load starting location
seeds = np.loadtxt(seedfile)
if len(seeds.shape) == 1:
    seeds = np.reshape(seeds, (1,10))
x0 = seeds[seedind,:]
if len(x0.shape) == 1 or x0.shape[0] == 1:
    p0 = emcee.utils.sample_ball(x0, len(x0)*[1e-8], nwalkers)
else:
    p0 = x0
    assert p0.shape[0] == nwalkers, "seed file should have nwalkers rows."

# Load observables
with open(obsfile, 'rb') as f:
    obs = pickle.load(f)

# Define log-likelihood function
def lnfake(x,n):
    return -0.5*np.dot(x,x)
def lnprob(x,n):
    (svec, dvec) = mdl(n, x, ref)
    P = pfh(x) + priors.gasplanet(svec, dvec, obs)
    if P > priors.whybother(): # usually whybother=-inf
        mrot = obs['m']
        Js,_ = tof4.tof4(svec, dvec, mrot, sskip=sskip)
        L = lfh(Js, obs)
        mL2 = ((dprofs.mass(svec, dvec) - obs['M'])/obs['dM'])**2
        L = np.sqrt(mL2 + L**2)
    else:
        L = 0.0 # note this does NOT mean a good score!
    return (-0.5*L**2) + P

# Open schwimmbad pool
if pooltype is 'serial':
    pool = schwimmbad.SerialPool()
elif pooltype is 'mpi':
    pool = schwimmbad.MPIPool()
else:
    raise ValueError("pooltype is 'serial' or 'mpi'")
if not pool.is_master():
    pool.wait()
    sys.exit(0)

# Define emcee ensemble sampler
if fakelike:
    sampler = emcee.EnsembleSampler(nwalkers, 10, lnfake, args=[toflevels], pool=pool)
else:
    sampler = emcee.EnsembleSampler(nwalkers, 10, lnprob, args=[toflevels], pool=pool)
sampler.a = ascale

# Run emcee ensemble sampler, in batches
nbatches = nsamples//batchsize
rtic = timer()
for batch in range(nbatches):
    cout("\nRunning batch {} of {}...".format(batch+1,nbatches))
    sys.stdout.flush()
    tic = timer()
    if batch is 0:
        sampler.run_mcmc(p0, batchsize)
    else:
        sampler.run_mcmc(None, batchsize)
    toc = timer()
    cout("done. ({:0.2g} hours, {:0.2g} sec per walker-sample)\n".format(
        (toc - tic)/3600, (toc - tic)/(batchsize + 1)))
    cout("Running mean acceptance ratio = {:0.2g}\n".format(
        np.mean(sampler.acceptance_fraction)))
    fname = "dump{}.pickle".format(batch%2)
    with open(fname,'wb') as f:
        sampler.pool = None
        lpfun = sampler.lnprobfn
        sampler.lnprobfn = None
        params = dict([('obs',obs), ('pfh',pfh), ('lfh',lfh), ('mdl',mdl),
                       ('toflevels',toflevels)])
        pickle.dump((sampler, params), f)
        sampler.lnprobfn = lpfun
        sampler.pool = pool
        cout("Sampler pickled in {}\n".format(fname))
    fname = "dump{}.npz".format(batch%2)
    np.savez_compressed(fname, C=sampler.chain)
    cout("Sampler state dumped in {}\n".format(fname))
cout("\nRun finished ({:0.2g} hours).\n".format((timer() - rtic)/3600))
cout("Acceptance ratio = {:0.2g}\n".format(np.mean(sampler.acceptance_fraction)))

# Close schwimmbad pool
pool.close()

# Save final sampler state
np.savez_compressed('final.npz',C=sampler.chain,L=sampler.lnprobability)
cout("Final state dumped in {}\n".format('final.npz'))
fname = "final.pickle"
with open(fname, 'wb') as f:
    sampler.pool = None
    sampler.lnprobfn = None
    params = dict([('obs',obs), ('pfh',pfh), ('lfh',lfh), ('mdl',mdl),
                   ('toflevels',toflevels)])
    pickle.dump((sampler, params), f)
    cout("Sampler pickled in {}\n\n".format(fname))
