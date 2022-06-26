###############################################################################
# Methods for working with the PPWD density parameterization
###############################################################################
import numpy as np

def ppwd_profile(N, x, rho0, zvec=None, forcemono=True):
    """Single polynomial plus two smoothed step functions.

    ppwd_profile(N, x, rho0) returns an N-point density profile, a (zvec,dvec)
    tuple, where density d is approximated by a single polynomial of normalized
    radius z=s/s0 overlain with two smoothed-step-functions. The default radii
    spacing is one of equal increments between z=1 and z=1/N. The parameter
    vector x and the surface density rho0 together fully define the density
    profile.

    The first three elements of x define the location, scale, and sharpness of
    the first (inner) step. Specify location 0<z_1<1 in normalized mean radius
    and specify scale in kg/m^3. The sharpness parameter is non-dimensional
    positive; a larger value results in a narrower step. Try ~100. The next
    three elements of x specify the same parameters for the second (outer)
    step. The remaining elements are the coefficients of the polynomial.
    Remember that a deg-n polynomial will be constrained by boundary conditions
    and defined with n-1 free coefficients, rather than n+1.

    You can BYO radii distribution with ppwd_profile(..., zvec), in which case
    zvec had better be a 1d vector and N will be ignored. The returned profile
    is defined on normalized radii, even if zvec wasn't.

    ppwd_profile(..., forcemono=True) forces the resulting density profile to
    be monotonically non increasing. This shouldn't be necessary but
    occasionally for high-resolution models and high-order polynomials there is
    an unwanted local min somewhere in the first few layers.
    """

    y = _fixpoly(x[6:], rho0)
    dprof = polynomial_profile(N, y, zvec, forcemono)
    ro0 = dprof[1][0]
    dprof = _add_density_jump(dprof, x[3], x[4], x[5]) # outer first
    dprof = _add_density_jump(dprof, x[0], x[1], x[2]) # inner last
    dvec = dprof[1] - dprof[1][0] + ro0
    zvec = dprof[0]
    return(zvec,dvec)

def polynomial_profile(N, x, zvec=None, forcemono=True):
    """Return density profile following single polynomial in normalized radius.

    polynomial_profile(N, x) returns an N-point density profile, i.e. a
    (zvec,dvec) tuple, where density d is approximated by a single polynomial
    of normalized radius z=s/s0, with coefficients x (in descending order). The
    default radii spacing is one of equal increments between z=1 and z=1/N.

    You can BYO radii distribution with polynomial_profile(N, x, zvec), in
    which case zvec had better be a 1d vector and N will be ignored. The
    returned profile is defined on normalized radii, even if zvec wasn't.

    polynomial_profile(..., forcemono=True) forces the resulting density
    profile to be monotonically non increasing. This is actually the default.
    It helps when a high degree polynomial tries hard to asymptote horizontally
    as x->1 and ends up with a tiny upward slope.
    """

    # Minimal input control
    assert np.isscalar(N)
    x = np.array(x)
    assert len(x.shape) == 1

    # This is a simple one
    if zvec is None:
        zvec = np.linspace(1, 1/N, N)
    else:
        zvec = zvec/zvec[0]
    dvec = np.polyval(x, zvec)
    if forcemono:
        dvec[0] = max(dvec[0], 0)
        for k in range(2,len(dvec)):
            dvec[k] = max(dvec[k], dvec[k-1])

    return (zvec, dvec)

def ppwd_prior(x,obs=None):
    """A prior on the PPWD sampling-space parameters.

    The first six parameters are sampling-space transformations of the
    smooth-step properties: location, scale, and sharpness. We are okay with
    uniform priors on the physical values, but we need to transform those to
    the correct prior on the sampling-space values. The rest of x just needs to
    be bounded.
    """
    lp = np.zeros(x.shape)
    outof = lambda a, S: a < min(S) or a > max(S)
    bexp = lambda y: min(np.exp(y), 1e307)
    logit = lambda y: np.log(y) - np.log(1 - y)
    
    if obs is None:
        class obs:
            rhomax = 3e4

    # inner z ~ U(lo,hi) -> y=logit(z) ~ (e^y)/(1 + e^y)^2*U(logit(lo,hi))
    support = (logit(0.05),logit(0.5))
    if outof(x[0], support):
        lp[0] = -np.inf
    else:
        lp[0] = x[0] - 2*np.log(1 + bexp(x[0]))

    # inner drho ~ U(0,rhomax) -> y=log(rhoc) ~ e^y*U(-inf,log(rhomax))
    support = (-20,np.log(obs.rhomax))
    if outof(x[1], support):
        lp[1] = -np.inf
    else:
        lp[1] = x[1]

    # inner sharpness ~ U(lo,hi) -> y=log(sharpness) ~ e^y*U(log(lo),log(hi))
    support = np.log((20,1001))
    if outof(x[2], support):
        lp[2] = -np.inf
    else:
        lp[2] = x[2]

    # outer z ~ U(lo,hi) -> y=logit(z) ~ (e^y)/(1 + e^y)^2*U(logit(lo,hi))
    support = (logit(0.5),logit(0.85))
    if outof(x[3], support):
        lp[3] = -np.inf
    else:
        lp[3] = x[3] - 2*np.log(1 + bexp(x[3]))

    # outer drho ~ U(0,rhomax) -> y=log(rhoc) ~ e^y*U(-inf,log(rhomax))
    support = (-20,np.log(obs.rhomax))
    if outof(x[4], support):
        lp[4] = -np.inf
    else:
        lp[4] = x[4]

    # outer sharpness ~ U(lo,hi) -> y=log(sharpness) ~ e^y*U(log(lo),log(hi))
    support = np.log((20,1001))
    if outof(x[5], support):
        lp[5] = -np.inf
    else:
        lp[5] = x[5]
    
    # polynomial coefficients, just a guess ~ U(-1e7,1e7)
    for k in range(6,len(x)):
        if outof(x[k], (-1e7,1e7)):
            lp[k] = -np.inf

    return sum(lp)

def ppwd_transform(x):
    """Transform mcmc sample vector to ppwd params."""
    y = np.concatenate((_fixcore(x[:3]), _fixcore(x[3:6]), x[6:]))
    return y

def _add_density_jump(dprof, z, scale, sharpness=100.0):
    """Add a localized density increase to existing profile.

    dprof = add_density_jump(dprof, z, scale) takes an existing N-point density
    profile, i.e., a (zvec,dvec) tuple, and adds a smoothed-step-function
    density localized in the neighborhood of normalized radius z. The scale
    parameter controls the overall height of the step (specified in real
    density units). In other words, scale*kg/m^3 will be added to the central
    density while nothing will be added to the surface density, with the bulk
    of the increase happening near z.

    add_density_jump(...,sharpness) allows you to control the step's
    "sharpness", the radial distance over which most of the density will be
    added. The sharpness parameter is dimensionless; experiment with different
    values to see the effect. Default 100.

    Algorithm: we use the inverse tangent function to approximate the step.
    """

    zvec, dvec = dprof
    x = zvec/zvec[0] # make sure radii are normalized before applying atan
    y = scale*(np.pi/2 + np.arctan(-sharpness*(x - z)))/np.pi
    dvec = dvec + y
    return (zvec,dvec)

def _fixpoly(x,rho0,z0=1.0,rho0p=None):
    """Return full polynomial coefficients from constrained ones.

    y = fixpoly(x,rho0) interprets the vector x as the n-1 coefficients of a
    degree-n polynomial that is constrained to have no linear term and pass
    through point (1,rho0). The returned vector will hold the n+1 coefficients
    of the full polynomial.

    y = fixpoly(x,rho0,z0) constrains that polynomial to pass through point 
    (z0,rho0) instead.

    y = fixpoly(x,rho0,z0,rho0p) additionally constrains the derivative to pass
    through (z0,rho0p). Note that this means a vector x of length n defines a
    polynomial of degree n+2.
    """

    if rho0p is None:
        return np.hstack((x,0,rho0-np.polyval(np.hstack((x,0,0)),1)))
    else:
        y = np.hstack((x,0,0,0))
        a2 = (rho0p - np.polyval(np.polyder(y),z0))/(2*z0)
        a0 = rho0 - np.polyval(np.hstack((x,a2,0,0)),z0)
        return np.hstack((x,a2,0,a0))

def _fixcore(x):
    """Transform core parameters from sampling space to physical space."""

    x = np.array(x)
    assert x.size == 3
    bexp = lambda y: min(np.exp(y), 1e307) # bounded exp to avoid inf/inf=nan
    expit = lambda y: bexp(y)/(1 + bexp(y))
    x[0] = expit(x[0]) # location, (-inf,inf)->(0,1)
    x[1] = bexp(x[1])  # scale (-inf,inf)->[0,inf)
    x[2] = bexp(x[2])  # sharpness (-inf,inf)->[0,inf)
    return x

def _unfixcore(x):
    """Transform core parameters from physical space to sampling space."""
    x = np.array(x)
    assert x.size == 3
    logit = lambda y: np.log(y) - np.log(1 - y)
    x[0] = logit(x[0])
    x[1] = np.log(x[1])
    x[2] = np.log(x[2])
    return x

def _ppwdref_profile(N, x, ref, zvec=None, forcemono=False):
    """(OBSOLETED) Referenced polynomial plus two smoothed step functions.

    ppwdref_profile(N, x, ref) returns an N-point density profile, a
    (zvec,dvec) tuple, where density d is approximated by a single polynomial
    of normalized radius z=s/s0 overlain with two smoothed-step-functions. The
    default radii spacing is one of equal increments between z=1 and z=1/N. The
    parameter vector x and the density array ref together fully define
    the density profile.

    The ref array is a reference density profile used as an extended boundary
    condition. The reference profile usually extends from z=1 down a short way
    (say to z=0.9) and the returned output dvec will match ref in that interval
    and also match the slope of ref at the lower end. It is assumed that ref is
    a proper profile (ref[1,1]==1, diff(ref[:,1])<0 and diff(ref[:,2])>0).

    The first three elements of x define the location, scale, and sharpness of
    the first (inner) step. Specify location 0<z_1<1 in normalized mean radius
    and specify scale in kg/m^3. The sharpness parameter is non-dimensional
    positive; a larger value results in a narrower step. Try ~100. The next
    three elements of x specify the same parameters for the second (outer)
    step. The remaining elements are the coefficients of the polynomial.
    Remember that a deg-n polynomial will be constrained by boundary conditions
    and defined with n-2 coefficients rather than n+1.

    You can BYO radii distribution with ppwdref_profile(..., zvec), in which
    case zvec had better be a 1d vector and N will be ignored. It doesn't
    matter if the input radii are normalized or not.

    ppwdref_profile(..., forcemono=True) forces the resulting density profile
    to be monotonically non increasing. This shouldn't be necessary but
    occasionally for high-resolution models and high-order polynomials there is
    an unwanted min somewhere in the first few layers.
    """

    if zvec is None:
        zvec = np.linspace(1, 1/N, N)
    zref = ref[:,0]
    dref = ref[:,1]
    rho_ref = dref[-1]
    rho_prime_ref = (dref[-2] - dref[-1])/(zref[-2] - zref[-1])
    y = fixpoly(x[6:], rho_ref, zref[-1], rho_prime_ref)
    dprof = polynomial_profile(N, y, zvec, forcemono)
    refind = dprof[0] >= zref[-1]
    dprof[1][refind] = np.interp(dprof[0][refind],zref[-1::-1],dref[-1::-1])
    ro0 = dprof[1][0]
    dprof = add_density_jump(dprof, x[3], x[4], x[5]) # outer first
    dprof = add_density_jump(dprof, x[0], x[1], x[2]) # inner last
    dvec = dprof[1] - dprof[1][0] + ro0
    return(zvec,dvec)
