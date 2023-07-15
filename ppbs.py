###############################################################################
# Methods for working with the PPBS barotrope parameterization
###############################################################################
from collections import namedtuple
import numpy as np
from TOFPlanet import TOFPlanet

_ppbs_supports = namedtuple("_ppbs_supports","K1 n1 K2 n2 K3 n3 z12 z23")
_def_sup = _ppbs_supports(
    n1=(0.02, 2.0),
    n2=(0.02, 2.0),
    n3=(0.02, 2.0),
    K1=(1e4, 8e5),
    K2=(1e4, 8e5),
    K3=(1e4, 8e5),
    z12=(0.45, 0.999),
    z23=(0.001, 0.8)
    )
_ben_sup = _ppbs_supports(
    n1=(0.02, 0.75),
    n2=(0.15, 1.3),
    n3=(0.2, 1.2),
    K1=(1e4, 8e5),
    K2=(8e4, 18e4),
    K3=(2.5e4, 1.2e5),
    z12=(0.45, 0.999),
    z23=(0.001, 0.8)
    )
_mono_y_seed = np.array([589361.3,
                         1.2993193778,
                         589361.3,
                         1.2993193778,
                         589361.3,
                         1.2993193778,
                         0.6,
                         0.3])
_beno_y_seed = np.array([100000,
                         0.3525024609,
                         129322.860274182,
                         0.8648796926,
                         56143.3213601624,
                         0.8081149898,
                         0.90177,
                         0.550049])

def ppbs_planet(N, y, obs, toforder=4, xlevels=-1):
    """Create a TOFPlanet from ppbs parameters."""
    p = TOFPlanet(obs)
    p.opts['toforder'] = toforder
    p.opts['xlevels'] = xlevels

    # The ppbs parameterization
    K1, n1, K2, n2, K3, n3, r12, r23 = y

    # Make radius grid; snap to layer boundaries
    zvec = np.linspace(1, 1/N, N)
    tind = np.argmin(np.abs(zvec - r12))
    cind = np.argmin(np.abs(zvec - r23))
    zvec[tind] = r12
    zvec[cind] = r23

    # Initialize with quadratic density profile
    a = _quadratic_planet(obs.M, obs.a0)
    p.si = zvec*obs.a0
    p.rhoi = a*zvec**2 - a

    # The ppbs parameters define the radius-dependent barotrope
    def tripoly(P):
        rho = np.full_like(P,np.nan)
        rho[:tind] = (P[:tind]/K1)**(1/(1 + 1/n1))
        rho[tind:cind] = (P[tind:cind]/K2)**(1/(1 + 1/n2))
        rho[cind:] = (P[cind:]/K3)**(1/(1 + 1/n3))
        return rho
    p.set_barotrope(tripoly)

    p._params = y
    return p

def ppbs_prior_uniform(x):
    """A uniform prior on the PPBS sampling-space parameters.

    We use independently uniform prior on the polytrope parameters and
    transiton radii, transformed those to the equivalent priors on the
    sample-space values. Concretely, if
        X = log(Y - a) - log(b - Y)
    then
        Y = (a + b*e^X)/(1 + e^X)
    and
        Y ~ U(a,b)
    becomes
        X ~ (b - a)*(e^x)/(1 + e^x)^2.
    """
    lp = np.zeros_like(x)
    for k in range(len(lp)):
        a, b = _def_sup[k]
        lp[k] = x[k] - 2*np.log(1 + np.exp(x[k])) + np.log(b - a)
    return sum(lp)

def _transform(x, supports=_def_sup):
    """Transform mcmc sample space to ppbs params."""
    y = np.full_like(x, np.nan)
    for k in range(x.size):
        y[k] = _expit(x[k], *supports[k])
    return y

def _untransform(y, supports=_def_sup):
    """Transform ppbs params vector to sample space."""
    x = np.full_like(y, np.nan)
    for k in range(y.size):
        x[k] = _logit(y[k], *supports[k])
    return x

def _quadratic_planet(M, R, rho0=0):
    """Return a s.t. rho(r) = a*(r/R)^2 - a integrates to M."""
    return 5/2*rho0 - 15*M/8/np.pi/R**3

def _expit(x,a,b):
    return (a + b*np.exp(x))/(1 + np.exp(x))

def _logit(x,a,b):
    return np.log(x - a) - np.log(b - x)

####
if __name__ == '__main__':
    print("alo")

