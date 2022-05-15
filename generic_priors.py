################################################################################
# Prior LOG-probabilities to ADD to log-likelihood
################################################################################
import numpy as np

def _unlikely():
    return -np.inf

def domain(x,D=(0,1)):
    """Verify parameter values inside domain."""
    x = np.array(x, ndmin=1)
    lp = np.zeros(x.shape)
    outof = lambda a, S: a <= min(S) or a >= max(S)
    for k in range(len(x)):
        if outof(x[k], D):
            lp[k] = _unlikely()
    return sum(lp)

def nonphysical_density(dvec):
    """Verify finite nonnegative density."""
    if any(np.isinf(dvec)):
        return _unlikely()
    elif any(np.isnan(dvec)):
        return _unlikely()
    elif np.any(dvec < 0):
        return _unlikely()
    else:
        return 0.0

def nonincreasing_density(dvec, flipum=False):
    """Verify monotonically downward increasing density."""
    if flipum:
        dvec = np.flipud(dvec)
    if any(np.diff(dvec) < 0):
        return _unlikely()
    else:
        return 0.0

def max_density(dvec, obs, flipum=False):
    """Apply softened density upper bound."""
    if flipum:
        dvec = np.flipud(dvec)
    romax = dvec[-1]
    d = (max(romax, obs.rhomax) - obs.rhomax)/obs.rhobar
    return -0.5*d**2

def gasplanet(dvec, obs, flipum=False):
    return (nonphysical_density(dvec) +
            nonincreasing_density(dvec,flipum) +
            max_density(dvec,obs,flipum))