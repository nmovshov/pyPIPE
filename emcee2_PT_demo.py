# try running:
# time mpirun -n 5 python emcee_PT_demo.py
# remember: one mpi process will be master, n-1 will run lnprob
import sys
import time
import numpy as np
import schwimmbad
import emcee
from timeit import default_timer as timer
cout = sys.stdout.write

ndim = 2
nwalkers = 8
ntemps = 20
nsamples = 100
batchsize = 10

# Define the bi-modal posterior we wish to sample
mu1 = np.ones(2)
mu2 = -np.ones(2)
sigma1inv = np.diag([100.0,100.0])
sigma2inv = np.diag([100.0,100.0])
def logl(x):
    dx1 = x - mu1
    dx2 = x - mu2
    return np.logaddexp(-0.5*np.dot(dx1, np.dot(sigma1inv, dx1)),
                        -0.5*np.dot(dx2, np.dot(sigma2inv, dx2)))

def logp(x):
    time.sleep(0.01)
    return 0.0

# The usual schwimmbad blurb
pool = schwimmbad.MPIPool()
if not pool.is_master():
    pool.wait()
    sys.exit(0)

# Construct the PT sampler and initial position
sampler = emcee.PTSampler(ntemps, nwalkers, ndim, logl, logp, pool=pool)
x0 = [0,0]
p0 = emcee.utils.sample_ball(x0, len(x0)*[1e-8], ntemps*nwalkers)
p0 = np.reshape(p0, (ntemps,nwalkers,ndim))

# Run sampler in batches
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
cout("\nRun finished ({:0.2g} hours).\n".format((timer() - rtic)/3600))
cout("Acceptance ratio = {:0.2g}\n".format(np.mean(sampler.acceptance_fraction)))
pool.close()
