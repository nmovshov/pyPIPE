################################################################################
# Template module. Copy and modify before running.
# For MPI run with:
# mpirun -n N python filename.py
################################################################################
from __future__ import division
import os, sys
import numpy as np
import tof4
import dprofs
import losses
import priors
import observables
from pymultinest.solve import Solver
from timeit import default_timer as timer
from warnings import filterwarnings
#filterwarnings('ignore', category=RuntimeWarning, message='invalid value')
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cout = sys.stdout.write
cerr = sys.stderr.write

### MODIFY LINES IN THIS SECTION ###
obs = observables.Saturn
quartfile = r'Quartic_fit-Mank_suite0.txt'
rtran = float(sys.argv[1])/100
rcore = float(sys.argv[2])/100
nlivepts = 1000
toflevels = 2048
sskip = 32
fakelike = False
verbose = True
###

# Load reference top atmosphere quartic fit
quart = np.loadtxt(quartfile)
assert len(quart == 5)

# Define log-likelihood function
def lnfake(x):
    from time import sleep
    sleep(0.012)
    modes = 8*[0]
    sigs = 8*[1]
    x = np.array(x)
    return sum(-0.5*((x - modes)/sigs)**2)
def lnprob(x):
    n = toflevels
    x = np.append(x, rtran)
    x = np.append(x, rcore)
    (svec, dvec) = dprofs.pppp_quartop(n, x, 0.94, obs.s0, quart)
    mrot = obs.m
    Js,out = tof4.tof4(svec, dvec, mrot, sskip=sskip)
    svec = svec*obs.a0/(svec[0]*out.a0)
    dsqr = (losses.euclid_J24(Js, obs)**2 + losses.mass((svec,dvec), obs)**2)
    return (-0.5*dsqr)

# Define PyMultinest solver
class PPPPQuartop(Solver):
    def Prior(self, cube):
        if fakelike:
            return -1000 + 2000*cube
        else:
            return priors.MN_ppqt_fixed_rads(cube)
    def LogLikelihood(self, cube):
        if fakelike:
            return lnfake(cube)
        else:
            return lnprob(cube)
        pass
    pass

# Directory to store output
outdir = 'mnchains-{:2.0f}{:2.0f}'.format(rtran*100, rcore*100)
if rank == 0:
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    else:
        cout("WARNING: directory {} already exists, ".format(outdir))
        cout("files may be overwritten.\n")

# Run PyMultinest solver
rtic = timer()
solution = PPPPQuartop(n_dims=8,
                       n_live_points=nlivepts,
                       outputfiles_basename=outdir+'/',
                       verbose=verbose)
rtoc = timer() - rtic

# Final words
if rank == 0:
    with open(os.path.join(outdir,'resume.dat'),'r') as f:
        nevals = int(f.readlines()[1].split()[1])
    del f
    cout("\n\n")
    cout("Run finished ({:0.2g} hours).\n".format(rtoc/3600))
    cout("N = {} likelihood evaluations ".format(nevals))
    cout("({:0.2g} seconds per evaluation).\n".format(rtoc/nevals))
    cout("Effective sample size {}.\n".format(solution.samples.shape[0]))
