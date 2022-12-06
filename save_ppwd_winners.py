import sys
import os
import numpy as np
import ppwd
import ahelpers as ah

if len(sys.argv) == 1:
    print("Usage: python save_ppwd_winners.py filename outputname dof {seedfile}")
    sys.exit(0)

fname = sys.argv[1]
outname = sys.argv[2]
df = int(sys.argv[3])
seedfile = sys.argv[4] if len(sys.argv) > 4 else None

C,L = ah.load_chain(fname)
print()
print("C.shape = ", C.shape)
keepers = ah.winners(C,L,df,ppwd.ppwd_prior) # 2.5 sigma with n dof
print("keepers.shape = ", keepers.shape)
Z = C[keepers,-1,:]
if Z.ndim == 1:
    Z = Z.reshape((1,-1))
print("Z.shape = ", Z.shape)

if seedfile:
    S = np.loadtxt(seedfile)
    fakers = []
    for z in Z:
        if z in S:
            fakers.append(True)
        else:
            fakers.append(False)
    Z = Z[~np.array(fakers)]
    print(f"found and removed {sum(fakers)} fakers.")

print(f"appending {Z.shape} to {outname}.")
try:
    ah._append(outname, Z)
except:
    np.savetxt(outname, Z)
