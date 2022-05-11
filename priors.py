################################################################################
# Prior LOG-probabilities to ADD to log-likelihood
################################################################################
import numpy as np

def nonincreasing_density(zvec, dvec, obs, flipum=False):
    if flipum:
        dvec = np.flipud(dvec)
    if any(np.diff(dvec) < 0):
        return _unlikely()
    else:
        return 0.0

def _upperbound_density_contrast(zvec, dvec, obs, flipum=False):
    """Obsolete check on central density as multiple of bulk density."""
    if flipum:
        zvec = np.flipud(zvec)
        dvec = np.flipud(dvec)
    dro = np.hstack((dvec[0], np.diff(dvec)))
    m = sum(dro*zvec**3)
    robar = m/zvec[0]**3
    dmax = dvec[-1]/robar
    romax = dmax*obs.rhobar
    d = (max(romax, obs.rhomax) - obs.rhomax)/obs.rhobar
    return -0.5*d**2

def rhomax(zvec, dvec, obs, flipum=False):
    if flipum:
        zvec = np.flipud(zvec)
        dvec = np.flipud(dvec)
    romax = dvec[-1]
    d = (max(romax, obs.rhomax) - obs.rhomax)/obs.rhobar
    return -0.5*d**2

def nonphysical_density(dvec):
    if any(np.isinf(dvec)):
        return _unlikely()
    elif any(np.isnan(dvec)):
        return _unlikely()
    elif np.any(dvec < 0):
        return _unlikely()
    else:
        return 0.0

def gasplanet(zvec, dvec, obs, flipum=False):
    if nonphysical_density(dvec) == 0.0:
        return nonincreasing_density(zvec, dvec, obs, flipum) + \
               rhomax(zvec, dvec, obs, flipum)
    else:
        return _unlikely()

def background_eos(zvec, dvec, obs, ss, SS, bgisen_list, flipum=False):
    pvec = _hydrostatic_pressure(zvec, dvec, obs, ss, SS, flipum)

    return pvec

def polynomial(x):
    # x = polynomial coefficients in descending order
    return 0

def fake_prior(x):
    x = np.array(x)
    lp = np.zeros(x.shape)
    outof = lambda a, S: a <= min(S) or a >= max(S)
    for k in range(len(x)):
        if outof(x[k], (0,1)):
            lp[k] = _unlikely()
    return sum(lp)

def whybother():
    return -np.inf

def _hydrostatic_pressure(zvec, dvec, obs, ss, SS, flipum=False):
    if flipum:
        zvec = np.flipud(zvec)
        dvec = np.flipud(dvec)

    # let's get the potential
    from tof4 import Upu
    N = len(ss[0]) - 1
    s0 = ss[0][N]; s2 = ss[1][N]; s4 = ss[2][N]; s6 = ss[3][N]; s8 = ss[4][N]
    aos = 1 + s0 - (1/2)*s2 + (3/8)*s4 - (5/16)*s6 + (35/128)*s8
    U = Upu(ss, SS, obs.m)
    U = np.flipud(U)  # WARNING: U was bottom up b/c of ss and SS
    U = U*obs.G*obs.M*aos/obs.a0*zvec**2

    # need density in real units now
    svec = zvec*obs.a0/aos
    drho = np.hstack((dvec[0],np.diff(dvec)))
    m = (4*np.pi/3)*sum(drho*svec**3)
    rvec = dvec*obs.M/m

    # integrate hydrostatic equilibrium top down
    rho = 0.5*(rvec[:-1] + rvec[1:])
    P = np.zeros(dvec.shape)
    P[0] = obs.P0
    P[1:] = P[0] + np.cumsum(-rho*np.diff(U))

    return P

def _unlikely():
    return -np.inf

def _test():
    x = [1,2,3]
    print(polynomial(x))

if __name__ == '__main__':
    _test()
