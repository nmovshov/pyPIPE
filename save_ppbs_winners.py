#------------------------------------------------------------------------------
# Script to select winners from a ppbs-sample emcee chain.
#------------------------------------------------------------------------------
import sys, os
import numpy as np
import argparse

# pyPIPE modules
import observables
import tof4, tof7
import TOFPlanet
import ahelpers as ah

import ppbs
the_mdl = ppbs.ppbs_planet
the_transform = ppbs._transform
the_prior = ppbs.ppbs_prior_uniform

def _main(args):
    # Load the raw emcee chain
    C,L = ah.load_chain(args.chainfile)
    print()
    print("C.shape = ", C.shape)

    # Candidates have loglike value minus prior equiv to 2.5 sigma with n dof
    keepers = ah.ppbs_winners(C,L,args.dof,the_prior,args.fix_rot)
    print("keepers.shape = ", keepers.shape)
    Z = C[keepers,-1,:]
    if Z.ndim == 1:
        Z = Z.reshape((1,-1))
    print("Z.shape = ", Z.shape)

    # Some of them could be fakers from original seed
    if args.seedfile and Z.size:
        S = np.loadtxt(args.seedfile)
        fakers = []
        for z in Z:
            if z in S:
                fakers.append(True)
            else:
                fakers.append(False)
        Z = Z[~np.array(fakers)]
        print(f"found and removed {sum(fakers)} fakers.")

    # In strict mode we cook the candidates and require independent scores
    if args.J_strict and Z.size:
        # Load planet observables
        try:
            obs = getattr(observables, args.observables)
        except:
            print("Could not find planet observables; check spelling?")
            sys.exit(1)

        colabs = []
        ind = [o//2 for o in range(2,args.J_strict+1,2)]
        for z in Z:
            p = ah.cook_ppbs_planet(z,obs,the_mdl,the_transform,**args.__dict__)
            E = np.abs(p.Js[ind] - obs.Js[ind])/obs.dJs[ind]
            if np.any(E > args.J_thresh):
                colabs.append(True)
            else:
                colabs.append(False)
        Z = Z[~np.array(colabs)]
        print(f"found and removed {sum(colabs)} collaborators.")

    # abyu
    print(f"appending {Z.shape} to {args.outname}.")
    try:
        ah._append(args.outname, Z)
    except:
        np.savetxt(args.outname, Z)

def _PCL():
    # Return struct with command line arguments as fields.
    parser = argparse.ArgumentParser(
        description="Select winners from tip of emcee chain.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('chainfile', help="npz file with emcee chain")
    parser.add_argument('outname', help="where to append winners to")
    parser.add_argument('dof', type=int,
        help="how many degrees of freedom in log-like")

    parser.add_argument('-s','--seedfile',
        help="Seed file used by driver")

    parser.add_argument('--J-strict', type=int, default=0,
        help="passing positive J-strict will require winners to " +
        "satisfy individual J values in addition to chi2 threshold")

    parser.add_argument('--J-thresh', type=float, default=2.5,
        help="a sigma-threshold for strict J conditions")

    parser.add_argument('-o','--observables',
            help="Observables struct (usually planet name, " +
                 "see observables.py for options")

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

    tofgroup.add_argument('--xlevels', type=int, default=128,
        help="Skip-n-spline levels")

    args = parser.parse_args()
    return args

if __name__ == '__main__':
    clargs = _PCL()
    _main(clargs)
