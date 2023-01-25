#------------------------------------------------------------------------------
# Drive emcee's ensemble sampler with our ppwd model. In this model the
# planet-polynomial is missing a linear term and constrained with a rho0
# boundary condition at z=1, and thus has n-1 free parameters for a polynomial
# of degree n. In addition, two arc-tan density jumps are modeled with three
# parameters each. So, for example, a degree-8 ppwd has 13 free parameters.
#
# Run this driver with:
# python drive-ppwd.py seedfile nsteps obs
#
# The target planet is determined by the third parameter: obs is the name of a
# class in the observables module, e.g., Neptune.
#
# There are many command-line arguments that control both modeling and sampling
# behavior. Run python drive-ppwd.py --help for the full list.
#
# Most tof options can be left at default. Use the _PCL() function to
# introduce more command-line arguments if needed.
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

# pyPIPE modules
import losses
import generic_priors
import observables
import tof4, tof7
import TOFPlanet
import ahelpers as ah

### PROJECT-SPECIFIC COMPONENTS ###
import ppwd
the_prior = ppwd.ppwd_prior
the_mdl = ppwd.ppwd_profile
the_transform = ppwd.ppwd_transform
###

### GENERIC PIPE COMPONENTS ###
# The rest of the driver defines the required components of a sampler and takes
# care of initializing the ensemble sampler and the parallel pools and
# batching/checkpointing system.
###

def _lnprob(x,obs,args):
    """Generate, relax, and evaluate a model TOFPlanet."""
    # Parameter vector x may or may not include (standardized) rotation period
    if args.fix_rot:
        Prot = obs.P
        xx = x
    else:
        Prot = x[0]*obs.dP/2 + obs.P
        xx = x[1:]

    # Evaluate prior on sample-space parameters
    P = (the_prior(xx, obs) + generic_priors.rotation_prior(Prot, obs))
    if P == -np.inf:
        return P

    # Transform from sample space to model space and create the density profile
    y = the_transform(xx, obs)
    svec, dvec = the_mdl(args.toflevels, y, obs.rho0)
    svec = svec*obs.s0

    # Evaluate prior on pre-solved planet
    if not args.fakelike:
        P = P + generic_priors.gasplanet(dvec, obs)

    # If model not pre-rejected, relax to HE and evaluate
    dsqr = 0
    if (P > generic_priors._unlikely()) and (not args.fakelike):
        tp = TOFPlanet.TOFPlanet(obs,**args.__dict__)
        tp.si = svec
        tp.rhoi = dvec
        tp.period = Prot
        if args.preserve_period:
            tp.relax_to_rotation()
        else:
            tp.mrot = (2*np.pi/tp.period)**2*tp.s0**3/tp.GM
            tp.relax_to_HE(fixmass=args.fix_mass,moi=args.with_moi)
        svec, dvec = tp.si, tp.rhoi
        Js = tp.Js

        jflag = args.Jays[args.Jays > 0]
        dsqr = (losses.mass((svec,dvec),obs)**2 +
                losses.rho0((svec,dvec),obs)**2 +
                losses.euclid_Jnm(Js,obs,jflag)**2)
        if args.with_moi:
            dsqr += losses.NMoI(tp.NMoI, obs)**2
        if args.with_k2:
            dsqr += losses.k2((svec, dvec), obs)**2
    return -0.5*(dsqr/args.temperature) + P

def _read_seeds(args, outdir, obs):
    """Generic seed reader.
    
    In pyPIPY drivers, args.seedfile is a usually csv file with a line for each
    seed. If there is a single seed the ensemble sampler will be initiated with
    a random perturbation of it for each walker. If there are as many seeds as
    walker, each walker is assigned a seed. If there are more seeds than
    walkers, the random subset will be chosen. If there are fewer seeds than
    walkers, but more than one, this is an error.

    However, we first look for a reseeds.txt file in the output directory as a
    simple hack to restart from previous run.
    """
    try:
        seedfile = os.path.join(outdir, 'reseeds.txt')
        if os.path.isfile(seedfile):
            p0 = np.loadtxt(seedfile, delimiter=',')
            if  len(p0.shape) == 1:
                p0 = np.reshape(p0, (1,-1))
        else:
            seedfile = args.seedfile
            try:
                p0 = np.loadtxt(seedfile, delimiter=',')
            except:
                p0 = np.loadtxt(seedfile, delimiter=' ')
            if len(p0.shape) == 1: # else it's on the user
                p0 = np.reshape(p0, (1,-1))
        
        cout("\nSeeds loaded from " + seedfile + ".\n")
    except:
        print("Failed to load seed(s); check command line.")
        raise

    np.savetxt(os.path.join(outdir,'seed.txt'), p0)
    return p0

def _main(spool,args):
    # Load and customize planet observables
    try:
        obs = getattr(observables, args.observables)(
                J2_sig=args.J2_error,
                J4_sig=args.J4_error,
                M_sig=args.M_error,
                dk2=args.k2_error)
    except:
        raise ValueError(
                "Could not determine target planet; check command line.")
    if args.no_spin: # Sometimes (rarely) we want to stop rotation
        obs.P = np.inf
        obs.m = 0.0
        obs.a0 = obs.s0

    # Make a directory to store output
    outdir = '{}_{}_run'.format(args.prefix,obs.pname)
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    else:
        cout("\nWARNING: directory {} already exists, ".format(outdir))
        cout("files may be overwritten.\n")

    # Load starting positions
    if args.restart:
        p0 = None
        print("Restarting from {}/state.h5.".format(outdir))
    else:
        p0 = _read_seeds(args, outdir, obs)
    
    # Determine dimensions and seed walkers
    if not args.restart:
        ndims =  p0.shape[1]
        if args.nwalkers is None:
            nwalkers = 2*ndims
        else:
            nwalkers = args.nwalkers
        if p0.shape[0] == 1: # starting fresh from single seed
            p0 = np.random.normal(p0, 1e-8, size=(nwalkers,ndims))
        elif p0.shape[0] == nwalkers: # maybe a restart from reseeds.txt
            p0 = p0
        elif p0.shape[0] > nwalkers:
            p0 = np.random.permutation(p0)[:nwalkers,:]
        else:
            print("ERROR: seed file does not have enough rows.")
            sys.exit(0)

    # Define the emcee ensemble sampler
    if args.serialize:
        backend = emcee.backends.HDFBackend(os.path.join(outdir,'state.h5'))
    else:
        backend = emcee.backends.Backend()
    if args.restart:
        C = backend.get_chain()
        nsteps_completed = C.shape[0]
        nwalkers = C.shape[1]
        ndims = C.shape[2]
        args.nsteps = max(args.nsteps - nsteps_completed, 0)
        print(f"{nsteps_completed} completed steps found in {outdir}/state.h5")
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

    sampler = emcee.EnsembleSampler(nwalkers, ndims, _lnprob, args=[obs,args],
            pool=pool, backend=backend, moves=moves)
    print("Using {} walkers on {} dimensions.".format(nwalkers, ndims))
    print("Temperature T={}.".format(args.temperature))

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
            cout("done. ({:0.2g} hours, {:0.2g} sec per walker-step)\n".format(
                (toc - tic)/3600, (toc - tic)/(batchsize)))
            cout("Running mean acceptance ratio = {:0.2g}\n".format(
                np.mean(sampler.acceptance_fraction)))
        fname = "dump{}.npz".format(batch%2)
        fname = os.path.join(outdir, fname)
        np.savez_compressed(fname, C=sampler.chain, L=sampler.lnprobability)
        if args.verbosity >= 1:
            cout("Chains dumped in {}\n".format(fname))

    cout("\nRun finished ({:0.2g} hours).\n".format((timer() - rtic)/3600))

    if args.verbosity >= 2:
        cout("Mean acceptance ratio = {:0.2g}\n".format(
            np.mean(sampler.acceptance_fraction)))
        cout("Mean autocorrelation time: {:0.2g} steps\n".format(
            np.mean(sampler.get_autocorr_time(quiet=True))))
        cout("Mean Gelman-Rubin score: {:0.4g}\n".format(
            np.mean(ws.ah._grubin(sampler.chain))))

    # Save final sampler state
    fname = os.path.join(outdir, 'final.npz')
    np.savez_compressed(fname, C=sampler.chain, L=sampler.lnprobability)
    cout("Final state dumped in {}\n".format(fname))

    # ABYU
    return

def _PCL():
    # Return struct with command line arguments as fields.

    parser = argparse.ArgumentParser(
        description="PPWD pyPIPE driver with emcee.",
        epilog="NOTE: For MPI run with mpirun -n N python driver.py"+
            " seed steps --mpi...",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('seedfile', help="Seed(s) to start sampling from")

    parser.add_argument('nsteps', type=int,
        help="Number of steps to be taken by each walker")

    parser.add_argument('observables',
            help="Specify observables struct (see observables.py for options")

    parser.add_argument('-v', '--verbosity', type=int, default=1,
        help="Control runtime message verbosity")

    parser.add_argument('-r', '--restart', action='store_true',
        help="Restart using emcee.Backend")

    parser.add_argument('-s', '--serialize', action='store_true',
        help="Serialize using emcee.Backend")

    parser.add_argument('-p', '--prefix', default='',
        help="Base name for output directory")

    likegroup = parser.add_argument_group('Likelihood options')

    likegroup.add_argument('-T', '--temperature', type=float, default=1.0,
        help="Simulated annealing temperature.")

    likegroup.add_argument('-j', '--Jays', type=int, nargs='*', default=[2,4],
        help="J-coefficients to include in likelihood; 0 to ignore gravity")

    likegroup.add_argument('--J2-error', type=float, default=None,
        help="obs.dJ2 will be this multiplied by obs.J2")

    likegroup.add_argument('--J4-error', type=float, default=None,
        help="obs.dJ4 will be this multiplied by obs.J4")

    likegroup.add_argument('--M-error', type=float, default=None,
        help="obs.dM will be this multiplied by obs.M")

    likegroup.add_argument('-I','--with-moi', action='store_true',
        help="Use NMoI in likelihood evaluation")

    likegroup.add_argument('--with-k2', action='store_true',
        help="Use tidal Love number K2 in likelihood evaluation")

    likegroup.add_argument('--k2-error', type=float, default=None,
        help="Overrides obs.dk2")

    likegroup.add_argument('-f', '--fakelike', action='store_true',
        help="Use fake (uniform) likelihood function (e.g. to test prior)")

    mdlgroup = parser.add_argument_group('Additional model options')

    mdlgroup.add_argument('--fix-mass', type=int, default=1,
        help="Normalize converged model mass to obs.M")

    mdlgroup.add_argument('--fix-rot', type=int, default=1,
        help="Don't sample rotation period (use obs.P instead)")

    mdlgroup.add_argument('--preserve-period', type=int, default=1,
        help="Iterate on rotation m until rotation period matched obs.P.")

    mdlgroup.add_argument('--no-spin', action='store_true',
        help="Make spherical planet (sets obs.P to inf)")

    tofgroup = parser.add_argument_group('TOF options',
        'Options controlling ToF gravity calculation')

    tofgroup.add_argument('--toforder', type=int, default=4, choices=[4,7],
        help="Theory of figures expansion order")

    tofgroup.add_argument('--toflevels', type=int, default=4096,
        help="Number of level surfaces used to discretize density profile")

    tofgroup.add_argument('--xlevels', type=int, default=256,
        help="Skip-n-spline levels")

    emceegroup = parser.add_argument_group('emcee options',
        'Additional options to control the emcee sampler')

    parser.add_argument('-w', '--nwalkers', type=int, default=None,
        help="Number of ensemble walkers; defaults to 2x the seed dimensions")

    emceegroup.add_argument('--ascale', type=float, default=2.0,
        help="Scale parameter for stretch moves")

    emceegroup.add_argument('--swalk', type=int, default=None,
        help="Number of helper-walkers for walk moves")

    emceegroup.add_argument('--moves', choices=['stretch','walk','de'],
        default='stretch',
        help="EXPERIMENTAL: ensemble move strategy")

    swimgroup = parser.add_mutually_exclusive_group()

    swimgroup.add_argument('--ncores', type=int, default=1,
        help="Use python multiprocessing")

    swimgroup.add_argument('--mpi', action='store_true',
        help="Run on multi cores (must call with mpirun -n > 1)")

    args = parser.parse_args()

    assert np.all(np.remainder(args.Jays,2) == 0), "Js must be even."
    assert np.all(np.array(args.Jays) >= 0), "Js must be nonnegative"
    assert args.temperature > 0, "Temperature must be positive."
    args.Jays = np.array(args.Jays)
    if args.restart:
        args.serialize = True
    if args.no_spin:
        args.Jays = np.array([0])
        args.fix_rot = 1

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
