#------------------------------------------------------------------------------
# Script to convert a sample to list of ppwd profiles.
#------------------------------------------------------------------------------
import sys, os
import numpy as np
import argparse
import schwimmbad
from timeit import default_timer as timer

# pyPIPE modules
import observables
import tof4, tof7
import ahelpers as ah
import ppwd

the_mdl = ppwd.ppwd_profile
the_transform = ppwd.ppwd_transform

class Planet:
    """Holds interior structure vectors."""
    mass   = 0 # reference mass
    radius = 0 # reference radius (equatorial!)
    period = 0 # reference rotation period
    P0     = 0 # reference pressure
    si     = 0 # vector of mean radii (top down, s0=si[0] is outer radius)
    rhoi   = 0 # vector of densities on si grid
    Js     = 0 # external gravity coefficients (returned by tof<n>)
    M      = 0 # calculated mass
    a0     = 0 # calculated equatorial radius
    s0     = 0 # surface mean radius (another name for si[0]
    mi     = 0 # cumulative mass below si
    ai     = 0 # equatorial radii on level surfaces
    rhobar = 0 # calculated mean density
    wrot   = 0 # rotation frequency, 2pi/period
    qrot   = 0 # rotation parameter wrot^2a0^3/GM
    mrot   = 0 # rotation parameter, wrot^2s0^3/GM
    aos    = 0 # calculated equatorial to mean radius ratio (from tof<n>)
    G = 6.67430e-11; # m^3 kg^-1 s^-2 (2018 NIST reference)

def cook_planet(x, obs, opts):
    """Create a planet object from sample-space parameters."""
    from time import sleep
    sleep(0.01)

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

    swimgroup = parser.add_argument_group('Parallel execution option')

    swimgroup.add_argument('--ncores', type=int, default=1,
        help="Use multiple cores on single node")

    args = parser.parse_args()
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
    tic = timer()
    for k in range(nsamp):
        if (args.ncores == 1) and (args.verbosity > 0):
            print(f"cooking planet {k+1} of {nsamp}...",end='')
        s = sample[k]
        p = cook_planet(s,obs,args)
        if (args.ncores == 1) and (args.verbosity > 0):
            print("done.")
    toc = timer()
    print(f"Cooking planets...done ({(toc-tic)/3600:0.2g} hours).")

if __name__ == "__main__":
    clargs = _PCL()

    # The schwimmbad stump
    with schwimmbad.choose_pool(mpi=False,processes=clargs.ncores) as pool:
        _main(pool,clargs)
