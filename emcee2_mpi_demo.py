# try running:
# time mpirun -n 5 python emcee_mpi_demo.py
# remember: one mpi process will be master, n-1 will run lnprob
import sys
import time
import numpy as np
import schwimmbad
import emcee
from timeit import default_timer as timer
cout = sys.stdout.write

ndim = 4
nwalkers = 8
nsamples = 100
batchsize = 10
p0 = [np.random.rand(ndim) for _ in range(nwalkers)]

def lnprob(x):
    time.sleep(0.1)
    return -0.5*np.dot(x,x)

pool = schwimmbad.MPIPool()
if not pool.is_master():
    pool.wait()
    sys.exit(0)

sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, pool=pool)

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
cout("\nRun finished ({:0.2g} hours).\n".format((timer() - rtic)/3600))
cout("Acceptance ratio = {:0.2g}\n".format(np.mean(sampler.acceptance_fraction)))
pool.close()
