#------------------------------------------------------------------------------
# Script to convert a sample to list of ppwd profiles.
#------------------------------------------------------------------------------
import sys, os
import numpy as np
import argparse
import schwimmbad

# pyPIPE modules
import observables
import tof4, tof7
import ahelpers as ah
import ppwd

the_mdl = ppwd.ppwd_profile
the_transform = ppwd.ppwd_transform

def _PCL():
    # Return struct with command line arguments as fields.
    parser = argparse.ArgumentParser(
        description="PPWD profiles from sample.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('samplefile', help="File with model params in rows")

    parser.add_argument('observables',
            help="Observables struct (usually planet name, " +
                 "see observables.py for options")

    parser.add_argument('--nsamples', type=int, default=None,
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
        except:
            sample = np.loadtxt(args.samplefile, delimiter=' ')
        if len(sample.shape) == 1:
            sample = np.reshape(sample, (1,-1))
        print(f"Found {sample.shape[0]} records in {args.samplefile}.")
    except:
        print("Failed to load sample; check command line.")
        raise

    # Load planet observables
    try:
        obs = getattr(observables, args.observables)
    except:
        print("Could not find planet observables; check spelling?")
        sys.exit(0)

if __name__ == "__main__":
    clargs = _PCL()

    # The schwimmbad stump
    with schwimmbad.choose_pool(mpi=False,processes=clargs.ncores) as pool:
        _main(pool,clargs)
