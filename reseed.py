#------------------------------------------------------------------------------
# To quickly restart a sampler, create a seed file from the last link of each
# walker and save back in the run dir. The driver can be set to load from
# reseeds.txt if found in the run dir.
# Run with:
#   python reseed.py rundir/final.npz
#------------------------------------------------------------------------------
import numpy as np
import sys
import os
import ahelpers as ah

if len(sys.argv) == 1:
    print("Usage: python reseed.py filename (usually rundir/final.npz)")
    sys.exit(0)

fname = sys.argv[1]
rundir,_ = os.path.split(fname)
seedfile = os.path.join(rundir, 'reseeds.txt')

C,_ = ah.load_chain(fname)
Z = np.squeeze(C[:,-1,:])

print()
print(f"Saving {Z.shape} seeds to {seedfile}.")
np.savetxt(seedfile, Z, delimiter=',')
