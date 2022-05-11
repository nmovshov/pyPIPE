#---------------------------------------------------------------------------------
#  Concentric Maclaurin Spheroids gravity calculator
#
# Author: Naor Movshovitz (nmovshov at gee mail dot com)
#---------------------------------------------------------------------------------
from __future__ import division
import sys
import numpy as np
import warnings
from numba import jit

def cms(zvec, dvec, qrot, dJtol=1e-6, maxiter=100, xlayers=64, J0s=None):
    """Return gravity coefficients of density profile in hydrostatic equilibrium.

    Parameters
    ----------
    zvec : ndarray, 1d, positive
        Equatorial radii of constant density layers.
    dvec : ndarray, 1d, positive
        Densities of layers. The density should be monotonically non-increasing
        with z, but this is not enforced.
    qrot : float, scalar, nonnegative
        Dimensionless rotation parameter, w^2*a^3/GM.
    xlayers : integer, nonnegative
        number of layers whose shape is explicitly calculated (rest will spline).
    J0s : struct
       J-like values representing initial state. This is not just for speeding up
       convergence. Mostly it's a mechanism to preserve state between calls.

    Outputs
    -------
    Js : 1-by-16 vector, real
        Even harmonic gravity coefficients J0 to J30. Typically only J2 to J10 are
        helpful. J0 is included as a sanity check and test of convergence.
    out : struct
        A structure holding other quantities calculated in the course of running
        cms. Including out.zetas and out.JLike that together define the converged
        hydrostatic shape.
    """

    # Minimal input control
    zvec = np.array(zvec)
    dvec = np.array(dvec)
    assert zvec.shape == dvec.shape
    assert qrot >= 0
    assert xlayers >= 0

    # Normalize and flip radii and densities (it's safe to normalize a normal)
    p = np.argsort(zvec)
    zvec = zvec[p]
    dvec = dvec[p]
    if zvec[0] == 0:
        zvec[0] = np.spacing(1)
    dro = np.hstack((dvec[-1], np.diff(np.flipud(dvec))))
    m = sum(dro*np.flipud(zvec)**3)
    robar = m/zvec[-1]**3
    zvec = np.flipud(zvec/zvec[-1])
    dvec = np.flipud(dvec/robar)

    ## Initialize local variables (in CMS notation)
    lambdas = zvec
    deltas = np.hstack((dvec[0], np.diff(dvec)))
    nlay = len(lambdas)
    nangles = 48
    kmax = 30

    # Initialize zetas as spherical
    zetas = np.ones((nlay, nangles))

    # Initialize J-like quantities for spherical planet
    if J0s is None:
        Jlike = _allocate_spherical_Js(nlay, kmax, lambdas, deltas)
    else:
        Jlike = J0s

    # Abscissas and weights for Gaussian quadrature
    mus, gws = _gauleg(0, 1, nangles)

    # Precompute Legendre polynomials for fixed colatitudes (for gauss quad)
    class Ps:
        pass
    Ps.Pnmu = np.zeros((kmax+1,nangles))
    Ps.Pnzero = np.zeros((kmax+1,))
    for k in range(kmax+1):
        Ps.Pnmu[k,:] = _Pn(k, mus)
        Ps.Pnzero[k] = _Pn(k, 0)

    # Precompute powers of ratios of lambdas
    lamratpow = np.nan*np.zeros((kmax+2,nlay,nlay))
    for ii in range(nlay):
        for jj in range(nlay):
            for kk in range(kmax+2):
                lamratpow[kk,ii,jj] = (lambdas[ii]/lambdas[jj])**(kk)

    # The loop (see Hubbard, 2012 and ./notes/CMS.pdf)
    Js = Jlike.Jn # J0=0 ensures at least one iteration
    sskip = max(nlay//xlayers, 1)
    eps = np.finfo(float).eps
    for iter in range(maxiter):
        # Update shape with current gravity
        new_zetas = _skipnspline_zetas(Jlike, Ps, lamratpow, qrot, zetas, lambdas, sskip)

        # Update gravity with current shape
        new_Jlike = _update_Js(lambdas, deltas, new_zetas, Ps, gws)
        new_Js = new_Jlike.Jn

        # Check for convergence of J0-J8 to terminate...
        dJs = np.abs((Js - new_Js)/(Js+eps))
        if np.all(dJs[0:4] < dJtol):
            break
        # ... or update to new values
        Jlike = new_Jlike
        Js = new_Js
        zetas = new_zetas

    # It's not necessarily terrible to reach maxiter, but we'd want to know
    if iter == (maxiter - 1):
        warnings.warn('Spheroid shapes may not be fully converged.')

    # Return
    Js = new_Js
    class out:
        pass
    out.dJs = dJs
    out.iter = iter
    out.zetas = new_zetas
    out.lambdas = lambdas
    out.deltas = deltas
    out.Jlike = new_Jlike
    out.mus = mus
    out.gws = gws
    out.Ps = Ps

    return (Js, out)

#--------------------------------
#  Helper functions
#--------------------------------
def _allocate_spherical_Js(nlay,nmom,lambdas,deltas):
    """Return a struct of J-like quantities initialized to a spherical planet."""

    class Js:
        pass
    Js.tilde = np.zeros((nlay,(nmom+1)))
    Js.tildeprime = np.zeros((nlay,(nmom+1)))
    Js.Jn = np.zeros((nmom//2+1,))
    den = sum(deltas*lambdas**3)
    Js.tilde[:,0] = -1*(deltas*lambdas**3)/den
    Js.tildeprime[:,0] = -1.5*(deltas*lambdas**3)/den
    Js.tildeprimeprime = 0.5*(deltas*lambdas**3)/den
    Js.Jn[0] = sum(Js.tilde[:,0])
    return Js

def _gauleg(x1, x2, n):
    """Return abscissas and weights for Gauss-Legendre n-point quadrature.

    x,w = GAULEG(x1,x2,n) returns the abscissas x and weights w that can be used
    to evaluate the definite integral, I, of a function well approximated by an
    (2n - 1) degree polynomial in the interval [x1,x2] using the Gauss-Legendre
    formula:
        I = sum(w.*f(x))

    Algorithm:
      This function is based on the C++ implementation of a routine with the
      same name in Numerical Recipes, 3rd Edition. But in several places I opt
      for readability over performance, on the assumption that this function is
      most likely to be called in a setup routine rather than in an inner-loop
      computation.

    Example:
      fun = np.sin
      x,w = gauleg(0, np.pi, 6)
      I_quad = scipy.integrate.quad(fun, 0, np.pi)
      I_gaussleg = sum(w*fun(x))

    Reference: William H. Press, Saul A. Teukolsky, William T. Vetterling, and
    Brian P. Flannery. 2007. Numerical Recipes 3rd Edition: The Art of Scientific
    Computing (3 ed.). Cambridge University Press, New York, NY, USA.
    """

    # Minimal assertions
    assert np.isfinite(x1), "x1 must be real and finite"
    assert np.isfinite(x2), "x2 must be real and finite"
    assert int(n) == n and n > 2, "n must be positive integer > 2"
    assert x2 > x1, "Interval must be positive"

    # Local variables
    tol = 1e-14
    m = int(np.ceil(n/2))
    xmid = (x1 + x2)/2
    dx = (x2 - x1)
    x = np.NaN*np.ones((n,))
    w = np.NaN*np.ones((n,))

    # Main loop
    for j in range(m):
        # Get j-th root of Legendre polynomial Pn, along with Pn' value there.
        z = np.cos(np.pi*(j + 0.75)/(n + 0.5)) # initial guess for j-th root
        while True:
            # Calculate Pn(z) and Pn-1(z) and Pn'(z)
            p = np.NaN*np.ones((n+1,))
            p[0] = 1
            p[1] = z
            for k in range(2,n+1):
                pkm1 = p[k-1]
                pkm2 = p[k-2]
                pk = (1/k)*((2*k - 1)*z*pkm1 - (k - 1)*pkm2)
                p[k] = pk
            pn = p[-1]
            pp = (n*p[-2] - n*z*p[-1])/(1 - z**2)

            # And now Newton's method (we are hopefully very near j-th root)
            oldz = z
            z = z - pn/pp
            if np.abs(z - oldz) < tol:
                break

        # Now use j-th root to get 2 abscissas and weights
        x[j]     = xmid - z*dx/2 # Scaled abscissa left of center
        x[n-1-j] = xmid + z*dx/2 # Scaled abscissa right of center
        w[j]     = dx/((1 - z**2)*pp**2)
        w[n-1-j] = w[j]

    # Verify and return
    assert np.all(np.isfinite(x))
    assert np.all(np.isfinite(w))
    return (x,w)

def _Pn(n, x):
    """Fast implementation of ordinary Legendre polynomials of low degree.

    y = Pn(n,x) returns the ordinary Legendre polynomial of degree n evaulated at
    x. For n <= 12 the polynomials are implemented explicitly resutling in faster
    calculation compared with the recursion formula. For n > 12 we fall back on
    scipy.special.eval_legendre(n,x).

    Note: in keeping with the premise of an optimized implementation this function
    performs no input checks at all. Use with care.
    """

    x = np.array(x, dtype=float)

    if n == 0:
        y = np.ones(x.shape)
    elif n == 1:
        y = x
    elif n == 2:
        y = 0.5*(3*x**2 - 1)
    elif n == 3:
        y = 0.5*(5*x**3 - 3*x)
    elif n == 4:
        y = (1/8)*(35*x**4 - 30*x**2 + 3)
    elif n == 5:
        y = (1/8)*(63*x**5 - 70*x**3 + 15*x)
    elif n == 6:
        y = (1/16)*(231*x**6 - 315*x**4 + 105*x**2 - 5)
    elif n == 7:
        y = (1/16)*(429*x**7 - 693*x**5 + 315*x**3 - 35*x)
    elif n == 8:
        y = (1/128)*(6435*x**8 - 12012*x**6 + 6930*x**4 - 1260*x**2 + 35)
    elif n == 9:
        y = (1/128)*(12155*x**9 - 25740*x**7 + 18018*x**5 - 4620*x**3 + 315*x)
    elif n == 10:
        y = (1/256)*(46189*x**10 - 109395*x**8 + 90090*x**6 - 30030*x**4 + 3465*x**2 - 63)
    elif n == 11:
        y = (1/256)*(88179*x**11 - 230945*x**9 + 218790*x**7 - 90090*x**5 + 15015*x**3 - 693*x)
    elif n == 12:
        y = (1/1024)*(676039*x**12 - 1939938*x**10 + 2078505*x**8 - 1021020*x**6 + 225225*x**4 - 18018*x**2 + 231)
    else:
        from scipy.special import eval_legendre
        y = eval_legendre(n,x)

    return y

def _zeta_j_of_alfa(j, alfa, Js, Ps, lamrats, qrot, oldzeta):
    """Call fzero on _eq52 to find equilibrium zeta."""
    from scipy.optimize import root_scalar
    Jt = Js.tilde
    Jtp = Js.tildeprime
    Jtpp = Js.tildeprimeprime
    P0 = Ps.Pnzero
    Pmu = Ps.Pnmu[:,alfa]
    def fun(x):
        return _eq52(x,j,Jt,Jtp,Jtpp,P0,Pmu,lamrats,qrot)
    if oldzeta == 1.0:
        x1 = oldzeta*0.99
    else:
        x1 = 1.0
    sol = root_scalar(fun, x0=oldzeta, x1=x1)
    return sol.root

@jit(nopython=True)
def _eq52(zja, jl, Jt, Jtp, Jtpp, P0, Pmu, lamrats, qrot):
    # locals
    nlay = lamrats.shape[1]
    kmax = len(P0)-1
    lamj3 = lamrats[3,jl,0]
    q = qrot
    pows = np.arange(kmax+2)
    zetpow = zja**pows
    zetipow = zja**-pows

    # first double sum
    y1 = 0
    for ii in range(jl,nlay):
        for kk in range(0,kmax,2):
            y1 = y1 + lamrats[kk,ii,jl]*Jt[ii,kk]*zetipow[kk]*Pmu[kk]

    # second double sum
    y2 = 0
    for ii in range(jl,nlay):
        for kk in range(0,kmax,2):
            y2 = y2 + lamrats[kk,ii,jl]*Jt[ii,kk]*P0[kk]

    # third double sum
    y3 = 0
    for ii in range(jl):
        for kk in range(0,kmax,2):
            y3 = y3 + lamrats[kk+1,jl,ii]*Jtp[ii,kk]*zetpow[kk+1]*Pmu[kk]

    # and forth double sum
    y4 = 0
    for ii in range(jl):
        for kk in range(0,kmax,2):
            y4 = y4 + lamrats[kk+1,jl,ii]*Jtp[ii,kk]*P0[kk]

    # a single sum
    y5 = 0
    for ii in range(jl):
        y5 = y5 + lamrats[3,jl,ii]*Jtpp[ii]*zetpow[3]

    # another single sum
    y6 = 0
    for ii in range(jl):
        y6 = y6 + lamrats[3,jl,ii]*Jtpp[ii]

    # and the rotation term
    y7 = -(1/3)*q*lamj3*zja**2*(1 - Pmu[2]) + (1/2)*q*lamj3

    # Now combine
    y = (1/zja)*(y1 + y3 + y5) - y2 - y4 - y6 + y7

    # And BYO
    return y

def _skipnspline_zetas(Js, Ps, lamrats, qrot, oldzetas, lambdas, sskip):
    """Update layer shapes using current value of Js."""
    # locals
    nlay = lamrats.shape[1]
    nangles = Ps.Pnmu.shape[1]
    ind = np.arange(0,nlay,sskip)
    Y = np.nan*np.zeros((len(ind),nangles))

    # skip...
    for j in range(len(ind)):
        for alfa in range(nangles):
            oldzeta = oldzetas[ind[j],alfa]
            Y[j,alfa] = _zeta_j_of_alfa(ind[j], alfa, Js, Ps, lamrats, qrot, oldzeta)

    #... and spline
    from scipy.interpolate import CubicSpline as spline
    newzetas = np.nan*np.zeros((nlay,nangles))
    for alfa in range(nangles):
        x = np.flipud(lambdas[ind])
        y = np.flipud(Y[:,alfa])
        cs = spline(x, y)
        newzetas[:,alfa] = cs(lambdas)

    assert np.all(np.isfinite(newzetas)), "interpolation error"
    return newzetas

def _update_Js(lambdas, deltas, zetas, Ps, gws):
    """Single-pass update of gravitational moments by Gaussian quad."""
    nlay = len(lambdas)
    kmax = len(Ps.Pnzero)-1

    # Precompute common denominator in eqs. (48)
    denom = 0
    for j in range(nlay):
        fun = zetas[j,:]**3
        I = sum(fun*gws) # the gauss quadrature formula
        denom = denom + I*deltas[j]*lambdas[j]**3

    # Do J tilde, eq. (48a)
    new_tilde = np.zeros((nlay,kmax+1))
    for ii in range(nlay):
        for kk in range(kmax+1):
            if (kk%2 > 0):
                continue
            fun = Ps.Pnmu[kk,:]*zetas[ii,:]**(kk+3)
            I = sum(gws*fun) # gauss quad formula
            new_tilde[ii,kk] = -(3/(kk + 3))*deltas[ii]*lambdas[ii]**3*I/denom

    # Do J tilde prime, eqs. (48b and 48c)
    new_tprime = np.zeros([nlay,kmax+1])
    for ii in range(nlay):
        for kk in range(kmax+1):
            if (kk%2 > 0):
                continue
            if kk == 2: # eq. (48c)
                fun = Ps.Pnmu[2,:]*np.log(zetas[ii,:])
                I = sum(gws*fun) # gauss quad formula
                new_tprime[ii,kk] = -3*deltas[ii]*lambdas[ii]**3*I/denom
            else:       # eq. (48b)
                fun = Ps.Pnmu[kk,:]*zetas[ii,:]**(2 - kk)
                I = sum(gws*fun) # gauss quad formula
                new_tprime[ii,kk] = -(3/(2 - kk))*deltas[ii]*lambdas[ii]**3*I/denom

    # Do J tilde double prime, eq. (48d)
    new_tpprime = np.zeros((nlay,))
    for ii in range(nlay):
        new_tpprime[ii] = 0.5*deltas[ii]*lambdas[ii]**3/denom

    # Return updated Js struct
    class newJs:
        pass
    newJs.tilde = new_tilde
    newJs.tildeprime = new_tprime
    newJs.tildeprimeprime = new_tpprime
    n = np.arange(0,kmax+1,2)
    newJs.Jn = np.zeros((len(n),))
    for k in range(len(n)):
        newJs.Jn[k] = sum(newJs.tilde[:,n[k]]*lambdas**n[k])

    return newJs

def _test():
    import time
    N = 12
    zvec = np.linspace(1, 1.0/N, N)
    dvec = np.linspace(1/N,2,N)
    qrot = 0.1
    tic = time.time()
    Js, out = cms(zvec,dvec,qrot,xlayers=12)
    toc = time.time()
    print(Js[0:3])
    print("Elapsed time {} seconds".format(toc-tic))

if __name__ == '__main__':
    _test()
    sys.exit(0)
