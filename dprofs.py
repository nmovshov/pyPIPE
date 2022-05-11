from __future__ import division
import sys
import numpy as np

def polynomial(N, x):
    """Return density profile following a single polynomial in normalized radius.

    polynomial(N, x) returns an N-point density profile, i.e. a (zvec,dvec) tuple,
    where density D is approximated by a single polynomial of normalized radoius
    z=s/s0, with coefficients x (in descending order).
    """

    # Minimal input control (some day)
    assert np.isscalar(N)
    x = np.array(x)
    assert len(x.shape) == 1

    # This is a simple one
    zvec = np.linspace(1, 1/N, N)
    dvec = np.polyval(x, zvec)
    return (zvec, dvec)

def tripquad_quartictop(N, x, z0, R_m, quart):
    """Return piecewise-quadratic profile stitched to a quartic atmosphere.

    tripquad_quartictop(N, x, z0, R_m, quart) returns an N-point density profile,
    i.e. a (zvec,dvec) tuple, where the density D is approximated by a
    piecewise-quadratic function in the normalized radius Z=r/R_m up to Z=z0
    and follows a quartic polynomial above z0. The array quart holds quartic
    coefficients in descending powers. Inside of z0, each quadratic segment is
    defined by its end points and a curvature coefficient. The segment breakpoints
    are at z1 (upper) and z2 (lower). The density at z0 is d10. The density at z1
    is d11 (the right-limit) and d21 the left-limit). The density at z2 is d22
    (the right-limit) and d32 (the left-limit). The density at z=0 is d33. The 11
    free parameters are ordered as follows:

        x = [a1, y10, y11, a2, y21, y22, a3, y32, y33, lz1, lz2]
    where
        a1: curvature of first segment (a1*x^2 + b1*x + c1)
        y10: d10
        y11: log(d11 - d10)
        a2: curvature of second segment (a2*x^2 + b2*x + c2)
        y21: log(d21 - d11)
        y22: log(d22 - d21)
        a3: curvature of third segment (a3*x^2 + b3*x + c3)
        y32: log(d32 - d22)
        y33: log(d33 - d32)
        lz1: log(z1) - log(1 - z1)
        lz2: log(z2) - log(1 - z2)

    NOTE: The densities must be specified in real units and the resulting model
    should have approximately the correct mass without need for renormalization.
    This is because the reference part of the planet (above Z0) cannot change.
    """

    # Minimal input control (some day)
    assert np.isscalar(N)
    assert np.isscalar(z0)
    assert np.isscalar(R_m)
    x = np.array(x)
    quart = np.array(quart)
    assert len(x.shape) == 1
    assert len(quart.shape) == 1

    # Local variables
    zvec = np.linspace(1, 1/N, N)
    svec = zvec*R_m

    # Interpreting the parameters
    a1 = x[0]; y10 = x[1]; y11 = x[2]
    a2 = x[3]; y21 = x[4]; y22 = x[5]
    a3 = x[6]; y32 = x[7]; y33 = x[8]
    lz1 = x[9]; lz2 = x[10]

    bexp = lambda x: min(np.exp(x), 1e307) # bounded exp to avoid inf/inf=nan
    d10 = y10
    d11 = d10 + bexp(y11)
    d21 = d11 + bexp(y21)
    d22 = d21 + bexp(y22)
    d32 = d22 + bexp(y32)
    d33 = d32 + bexp(y33)
    r0 = z0
    rt = bexp(lz1)/(1 + bexp(lz1))
    rc = bexp(lz2)/(1 + bexp(lz2))

    # Upper envelope region
    b1 = (d11 - d10)/(rt - r0) - a1*(rt + r0)
    c1 = d10 - a1*r0**2 - b1*r0

    # Lower envelope region
    b2 = (d22 - d21)/(rc - rt) - a2*(rc + rt)
    c2 = d21 - a2*rt**2 - b2*rt

    # Core region
    b3 = (d32 - d33)/rc - a3*rc
    c3 = d33

    # Write the profile
    dvec = np.zeros(N)
    for k in range(N):
        if zvec[k] > r0:
            dvec[k] = d10*np.polyval(quart, zvec[k])
        elif zvec[k] > rt:
            dvec[k] = a1*zvec[k]**2 + b1*zvec[k] + c1
        elif zvec[k] > rc:
            dvec[k] = a2*zvec[k]**2 + b2*zvec[k] + c2
        else:
            dvec[k] = a3*zvec[k]**2 + b3*zvec[k] + c3

    # Return
    return (svec, dvec)

def mass(svec, dvec):
    """Return approximate mass integral."""
    from scipy.integrate import trapz
    return -4*np.pi*trapz(dvec*svec**2, x=svec)

def _test():
    print("alo world")
    zvec, dvec = polynomial(12, (-1, 0, 1))
    print(zvec)
    print(dvec)

if __name__ == '__main__':
    _test()
