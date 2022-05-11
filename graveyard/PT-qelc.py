################################################################################
# Template module. Copy and modify before running.
# For MPI use pooltype='mpi' and run with:
# mpirun -n N python filename.py
# Remember: one mpi process will be master, N-1 will run lnprob
################################################################################
from __future__ import division
import os, sys
import numpy as np
import emcee
import schwimmbad
import tof4
import dprofs
import losses
import priors
import observables
from timeit import default_timer as timer
from warnings import filterwarnings
#filterwarnings('ignore', category=RuntimeWarning, message='invalid value')
cout = sys.stdout.write
cerr = sys.stderr.write

### MODIFY LINES IN THIS SECTION ###
obs = observables.Saturn
seedfile = r'qelc-seed5520.txt'
quartfile = r'Quartic_fit-Mank_suite0.txt'
toflevels = 2048
batchsize = 1
nwalkers = 20
ntemps = 20
#seedind = range(ntemps*nwalkers)
seedind = 0
ascale = 2.0
sskip = 32
pooltype = 'mpi'
fakelike = False
###

# The command line
if len(sys.argv) < 2:
    print("Usage:\tpython PT-qelc nsamples")
    sys.exit(0)
nsamples = int(sys.argv[1])

# Load reference top atmosphere quartic fit
quart = np.loadtxt(quartfile)
assert len(quart == 5)

# Define the model and supporting functions
N = toflevels
logit = lambda x: np.log(x) - np.log(1 - x)
def qelc_to_tqq(x):
    return np.insert(x, 6, 0.0)
def mdl(x):
    y = qelc_to_tqq(x)
    lof = dprofs.tripquad_quartictop(N, y, 0.94, obs.s0, quart)
    return lof
ndims = 10
def pfh(x):
    return priors.tripquad_quartop(qelc_to_tqq(x))

# Load starting location
seeds = np.loadtxt(seedfile)
if len(seeds.shape) == 1:
    seeds = np.reshape(seeds, (1,ndims))
x0 = seeds[seedind,:]
if len(x0.shape) == 1 or x0.shape[0] == 1:
    p0 = emcee.utils.sample_ball(x0, len(x0)*[1e-8], nwalkers*ntemps)
else:
    p0 = x0
    assert p0.shape[0] == nwalkers*ntemps, "seed file does not have enough rows."
p0 = np.reshape(p0, (ntemps,nwalkers,ndims))

# Define log-likelihood function
def lnfake(x):
    from time import sleep
    sleep(0.01)
    return 0.0
def lnprob(x):
    (svec, dvec) = mdl(x)
    P = priors.gasplanet(svec, dvec, obs)
    if P > priors.whybother(): # usually whybother=-inf
        mrot = obs.m
        Js,out = tof4.tof4(svec, dvec, mrot, sskip=sskip)
        svec = svec*obs.a0/(svec[0]*out.a0)
        dsqr = losses.euclid_J24(Js, obs)**2 + losses.mass((svec,dvec), obs)**2
        return (-0.5*dsqr) + P # (P is most likely zero here)
    else:
        return P

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
    sampler = emcee.PTSampler(ntemps, nwalkers, ndims, lnfake, pfh, pool=pool)
else:
    sampler = emcee.PTSampler(ntemps, nwalkers, ndims, lnprob, pfh, pool=pool)
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
    cout("done. ({:0.2g} hours, {:0.2g} sec per temp-walker-sample)\n".format(
        (toc - tic)/3600, (toc - tic)/(batchsize + 1)))
    cout("Running mean acceptance ratio = {:0.2g}\n".format(
        np.mean(sampler.acceptance_fraction)))
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
