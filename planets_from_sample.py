#------------------------------------------------------------------------------
# Script to convert a sample to list of ppwd profiles.
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

def cook_planet(x, obs, opts):
    """Create a planet object from sample-space parameters."""
    y = the_transform(x, obs)
    svec, dvec = the_mdl(opts.toflevels, y, obs.rho0)
    p = TOFPlanet.TOFPlanet(obs)
    p.si = svec*obs.s0
    p.rhoi = dvec
    p.opts['toforder'] = opts.toforder
    p.opts['xlevels'] = opts.xlevels
    p.relax_to_HE(fixradius=True, moi=opts.with_moi, pressure=opts.with_pressure)
    return p

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

    mdlgroup.add_argument('--fix-rho0', type=int, default=1,
        help="Don't sample 1-bar density (use obs.rho0 instead)")

    mdlgroup.add_argument('--fix-mrot', type=int, default=1,
        help="Don't sample rotation parameter (use obs.m instead)")

    mdlgroup.add_argument('--no-spin', action='store_true',
        help="Make spherical planet (sets obs.m to zero)")

    tofgroup = parser.add_argument_group('TOF options',
        'Options controlling ToF gravity calculation')

    tofgroup.add_argument('--toforder', type=int, default=4, choices=[4,7],
        help="Theory of figures expansion order")

    tofgroup.add_argument('--toflevels', type=int, default=4096,
        help="Number of level surfaces used to discretize density profile")

    tofgroup.add_argument('--xlevels', type=int, default=256,
        help="Skip-n-spline levels")

    tofgroup.add_argument('-m','--with-moi', type=int, default=1,
        help="Include MoI calculation")

    tofgroup.add_argument('-P','--with-pressure', type=int, default=1,
        help="Include pressure calculation")

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
        p = cook_planet(s,obs,args)
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
