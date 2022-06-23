#------------------------------------------------------------------------------
# Example of driving a pyPIPE sample. Copy to project folder and customize, at
# a minimum, the model, likelihood, and prior, and almost surely the seed
# loader.
#
# HOW TO USE MULTIPLE COMPUTE CORES?
# Option 1: on a single, multi-core machine, include a --ncores=N argument to
#           the python call, e.g.
#    python driver.py seed nsteps --ncores=4
# Option 2: on a multi-node cluster, include --mpi argument to python AND run
#           with mpirun/mpiexec:
#    mpirun -n 60 python driver.py seed nsteps --mpi
#
# Remember: one mpi process will be master, N-1 will run lnprob, and emcee runs
# at most HALF the walkers in parallel.
#------------------------------------------------------------------------------
import sys, os
import numpy as np
import argparse
import emcee
import schwimmbad
from timeit import default_timer as timer
cout = sys.stdout.write
cerr = sys.stderr.write
os.environ['OMP_NUM_THREADS'] = '1' # recommended by dfm

# import .... # project-specific workspace and/or modules

def pfh(x):
    """CUSTOMIZE THIS PRIOR FOR THE PROJECT."""
    return -0.0

def mloss():
    """CUSTOMIZE A LOSS FUNCTION IF NEEDED."""
    return

def lnprob(x):
    """CUSTOMIZE THIS LOG-LIKELIHOOD FUNCTION FOR THE PROJECT."""
    # For example, this is a 2d gaussian but with some compute cost.
    from time import sleep
    sleep(0.01)
    return -0.5*np.dot(x,x) + pfh(x) # NOTE INCLUDE THE PRIOR!

def read_seeds(args, outdir):
    """CUSTOMIZE A SEED READER IF NEEDED."""
    try: # for this example the seed is a scalar read from CL and duplicated
        p0 = np.reshape(np.array(2*[float(args.seed)]), (1,2))
        print("Seed read from command line")
    except:
        print("Failed to load seed(s); check command line.")
        raise
    np.savetxt(os.path.join(outdir,'seed.txt'), p0)
    return p0


def _main(spool,args):
    # Make a directory to store output
    outdir = '{}'.format(args.prefix)
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    else:
        cout("WARNING: directory {} already exists, ".format(outdir))
        cout("files may be overwritten.\n")

    # Load starting positions (usually from text file and project-specific)
    if args.restart:
        p0 = None
        print("Restarting from {}/state.h5.".format(outdir))
    else:
        p0 = read_seeds(args, outdir)

    # Determine dimensionality and seed walkers
    if not args.restart:
        ndims = p0.shape[1]
        if args.nwalkers is None:
            nwalkers = 2*ndims
        else:
            nwalkers = args.nwalkers
        if p0.shape[0] == 1:
            p0 = np.random.normal(p0, 1e-8, size=(nwalkers,ndims))

    # Define the emcee ensemble sampler
    if args.fakelike:
        lfh = lnfake
    else:
        lfh = lnprob
    if args.serialize:
        backend = emcee.backends.HDFBackend(os.path.join(outdir,'state.h5'))
    else:
        backend = emcee.backends.Backend()
    if args.restart:
        C = backend.get_last_sample().coords
        nwalkers = C.shape[0]
        ndims = C.shape[1]
    else:
        backend.reset(nwalkers, ndims)

    if args.moves == 'stretch':
        moves = emcee.moves.StretchMove(a=args.ascale)
    elif args.moves == 'walk':
        moves = emcee.moves.WalkMove(s=args.swalk)
    elif args.moves == 'de':
        moves = [(emcee.moves.DESnookerMove(), 0.0),
                 (emcee.moves.DEMove(), 0.9),
                 (emcee.moves.DEMove(gamma0=1.0), 0.1)]
    else:
        print("ERROR: unknown move strategy.""")
        sys.exit(0)

    # If lnprob needs args beyond x use args=[] in constructor here
    sampler = emcee.EnsembleSampler(nwalkers, ndims, lfh,
            pool=pool, backend=backend, moves=moves)
    print("Using {} walkers on {} dimensions.".format(nwalkers, ndims))
    print("Using {} move strategy.".format(args.moves))

    # OPTIONAL TWEAKS TO SAMPLER HERE
    #sampler.a = ...

    # Run ensemble sampler in batches
    batchsize = max(args.nsteps//10, 1)
    nbatches = args.nsteps//batchsize
    rtic = timer()
    for batch in range(nbatches):
        cout("\nRunning batch {} of {}...".format(batch+1,nbatches))
        sys.stdout.flush()
        tic = timer()
        if batch == 0:
            sampler.run_mcmc(p0, batchsize, skip_initial_state_check=True)
        else:
            sampler.run_mcmc(None, batchsize, skip_initial_state_check=True)
        toc = timer()
        if args.verbosity >= 1:
            cout("done. ({:0.2g} hours, {:0.2g} sec per walker-sample)\n".format(
                (toc - tic)/3600, (toc - tic)/(batchsize)))
            cout("Running mean acceptance ratio = {:0.2g}\n".format(
                np.mean(sampler.acceptance_fraction)))
        fname = "dump{}.npz".format(batch%2)
        fname = os.path.join(outdir, fname)
        np.savez_compressed(fname, C=sampler.chain)
        if args.verbosity >= 1:
            cout("Chains dumped in {}\n".format(fname))

    cout("\nRun finished ({:0.2g} hours).\n".format((timer() - rtic)/3600))

    if args.verbosity >= 2:
        cout("Mean acceptance ratio = {:0.2g}\n".format(
            np.mean(sampler.acceptance_fraction)))
        cout("Mean autocorrelation time: {:0.2g} steps\n".format(
            np.mean(sampler.get_autocorr_time(quiet=True))))

    # Save final sampler state
    fname = os.path.join(outdir, 'final.npz')
    np.savez_compressed(fname, C=sampler.chain, L=sampler.lnprobability)
    cout("Final state dumped in {}\n".format(fname))

    # ABYU
    return

def lnfake(x):
    """Use this uniform prior to test the prior pfh."""
    return -0.0 + pfh(x)

def _PCL():
    # Return struct with command line arguments as fields.

    parser = argparse.ArgumentParser(
        description="Example pyPIPE emcee-driver: sample 2D standard normal.",
        epilog="NOTE: For MPI run with mpirun -n N python driver.py"+
            " seed steps ...",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('seed', help="Seed(s) to start sampling from. " +
        "This is usually the name of a text file. In this example it is " +
        "a simple number.")

    parser.add_argument('nsteps', type=int,
        help="Number of steps to be taken by each walker.")

    parser.add_argument('-v', '--verbosity', type=int, default=1,
        help="Control runtime message verbosity.")

    parser.add_argument('-w', '--nwalkers', type=int, default=None,
        help="Number of ensemble walkers; defaults to twice the seed dimensions.")

    parser.add_argument('-r', '--restart', action='store_true',
        help="Restart using emcee.Backend.")

    parser.add_argument('-s', '--serialize', action='store_true',
        help="Serialize using emcee.Backend")

    parser.add_argument('-p', '--prefix', default='run',
        help="Base name for output directory.")

    parser.add_argument('-f', '--fakelike', action='store_true',
        help="Use fake (uniform) likelihood function (e.g. to test prior).")

    emceegroup = parser.add_argument_group('emcee options',
        'Additional options to control the emcee sampler.')

    emceegroup.add_argument('--ascale', type=float, default=2.0,
        help="Scale parameter for stretch moves.")

    emceegroup.add_argument('--swalk', type=int, default=None,
        help="Number of helper-walkers for walk moves.")

    emceegroup.add_argument('--moves', choices=['stretch','walk','de'],
        default='stretch',
        help="EXPERIMENTAL: ensemble move strategy.")

    swimgroup = parser.add_mutually_exclusive_group()

    swimgroup.add_argument('--ncores', type=int, default=1,
            help="Use python multiprocessing.")

    swimgroup.add_argument('--mpi', action='store_true',
            help="Run on multi cores (must call with mpirun -n > 1).")

    args = parser.parse_args()

    if args.restart:
        args.serialize = True

    return args

if __name__ == "__main__":
    # Parse command line arguments
    clargs = _PCL()

    # The schwimmbad stump
    with schwimmbad.choose_pool(mpi=clargs.mpi,processes=clargs.ncores) as pool:
        if type(pool) is schwimmbad.multiprocessing.MultiPool:
            _main(pool,clargs)
        else:
            if not pool.is_master():
                pool.wait()
                sys.exit(0)
            _main(pool,clargs)
