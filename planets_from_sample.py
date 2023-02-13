#------------------------------------------------------------------------------
# Script to convert a ppwd-sample to list of TOFPlanets.
#------------------------------------------------------------------------------
import sys, os
import pickle
import numpy as np
import argparse
import schwimmbad
from timeit import default_timer as timer

# pyPIPE modules
import observables
import tof4, tof7
import TOFPlanet
import ahelpers as ah

import ppwd
the_mdl = ppwd.ppwd_profile
the_transform = ppwd.ppwd_transform

# def cook_planet(x, obs, opts):
#     """Create a planet object from sample-space parameters."""
#     if opts.fix_rot:
#         Prot = obs.P
#         x = x
#     else:
#         Prot = x[0]*obs.dP/2 + obs.P
#         x = x[1:]
#     y = the_transform(x, obs)
#     svec, dvec = the_mdl(opts.toflevels, y, obs.rho0)
#     tp = TOFPlanet.TOFPlanet(obs)
#     tp.si = svec*obs.s0
#     tp.rhoi = dvec
#     tp.period = Prot
#     tp.opts['toforder'] = opts.toforder
#     tp.opts['xlevels'] = opts.xlevels
#     if opts.no_spin:
#         tp.period = np.inf
#     if opts.preserve_period:
#         tp.relax_to_rotation(opts.fix_mass)
#     else:
#         tp.mrot = (2*np.pi/tp.period)**2*tp.s0**3/tp.GM
#     tp.relax_to_HE(fixmass=opts.fix_mass,moi=True,pressure=True)

#     if opts.with_k2:
#         tp.k2 = ah.lovek2(tp.si, tp.rhoi)
#     if not opts.savess:
#         tp.ss = None
#         tp.SS = None
#     return tp

def _PCL():
    # Return struct with command line arguments as fields.
    parser = argparse.ArgumentParser(
        description="PPWD profiles from sample.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('samplefile', help="File with model params in rows")

    parser.add_argument('observables',
            help="Observables struct (usually planet name, " +
                 "see observables.py for options")

    parser.add_argument('-n', '--nsamples', type=int, default=None,
        help="Number of samples to convert (leave blank to use entire sample)")

    parser.add_argument('-v', '--verbosity', type=int, default=1,
        help="Control runtime message verbosity")

    mdlgroup = parser.add_argument_group('Additional model options')

    mdlgroup.add_argument('--fix-mass', type=int, default=1,
        help="Normalize converged model mass to obs.M")

    mdlgroup.add_argument('--fix-rot', type=int, default=1,
        help="Don't sample rotation period (use obs.P instead)")

    mdlgroup.add_argument('--preserve-period', type=int, default=1,
        help="Iterate on rotation m until rotation period matched obs.P.")

    mdlgroup.add_argument('--no-spin', action='store_true',
        help="Make spherical planet (sets obs.P to inf)")

    mdlgroup.add_argument('--with-k2', action='store_true',
        help="Include tidal k2 calculation")

    tofgroup = parser.add_argument_group('TOF options',
        'Options controlling ToF gravity calculation')

    tofgroup.add_argument('--toforder', type=int, default=4, choices=[4,7],
        help="Theory of figures expansion order")

    tofgroup.add_argument('--toflevels', type=int, default=4096,
        help="Number of level surfaces used to discretize density profile")

    tofgroup.add_argument('--xlevels', type=int, default=256,
        help="Skip-n-spline levels")

    tofgroup.add_argument('--savess', action='store_true',
        help="Retain full shape information (triples file size)")

    swimgroup = parser.add_argument_group('Parallel execution option')

    swimgroup.add_argument('--ncores', type=int, default=1,
        help="Use multiple cores on single node")

    args = parser.parse_args()
    if args.ncores > 1:
        print("WARNING: multi-core support coming soon.")
        args.ncores = 1
    return args

def _main(spool,args):
    # Load sample file
    try:
        try:
            sample = np.loadtxt(args.samplefile, delimiter=',')
        except ValueError:
            sample = np.loadtxt(args.samplefile, delimiter=' ')
        if len(sample.shape) == 1:
            sample = np.reshape(sample, (1,-1))
        print(f"Found {sample.shape[0]} records in {args.samplefile}.")
    except:
        print("Failed to load sample; check command line.")
        sys.exit(1)

    # Load planet observables
    try:
        obs = getattr(observables, args.observables)
    except:
        print("Could not find planet observables; check spelling?")
        sys.exit(1)

    # Create a planet from each row
    nsamp = args.nsamples
    if nsamp is None:
        nsamp = sample.shape[0]
    print("Cooking planets...")
    planets = []
    tic = timer()
    for k in range(nsamp):
        if (args.ncores == 1) and (args.verbosity > 0):
            print(f"cooking planet {k+1} of {nsamp}...",end='')
        s = sample[k]
        p = ah.cook_planet(s,obs,the_mdl,the_transform,**args.__dict__)
        planets.append(p)
        if (args.ncores == 1) and (args.verbosity > 0):
            print("done.")
    toc = timer()
    print(f"Cooking planets...done ({(toc-tic)/3600:0.2g} hours).")

    # pickle planets
    outname = f"{args.samplefile[:-4]}_planets.pickle"
    with open(outname,'wb') as f:
        pickle.dump(planets,f)
    print(f"pickled {len(planets)} planets in {outname}.")

if __name__ == "__main__":
    clargs = _PCL()

    # The schwimmbad stump
    with schwimmbad.choose_pool(mpi=False,processes=clargs.ncores) as pool:
        _main(pool,clargs)
