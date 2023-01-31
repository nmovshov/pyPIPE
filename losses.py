###############################################################################
# Loss functions to pass to log-likelihood functions for sampling
###############################################################################
import numpy as np
import ahelpers as ah

def euclid_Jnm(Js, obs, ord=(2,4,6)):
    """Weighted euclidean distance in J space.

    This loss function is suitable for use when the J uncertainties are assumed
    to be uncorrelated. The dJ values in obs are interpreted as 1-sigma values
    for the purpose of weighing the distance. In other words, this is
    equivalent to the Mahalanobis distance with a diagonal covariance. A
    log-likelihood function proportional to a multivariate normal can simply
    use this negative one half times the square of this loss function.

    The tuple ord specifies which Js to use, in "math" not in python. So, e.g.,
    J2 is Js[1] and J4 is Js[2] etc. It is an error to specify orders that are
    missing from either Js or obs. For convenience, the default is ord=(2,4,6).
    """

    ind = [o//2 for o in ord]
    WD = np.array([(Js[i] - obs.Js[i])/obs.dJs[i] for i in ind])
    return np.sqrt(sum(WD**2))

def smooth_J_box(Js, obs, ord=(2,4,6)):
    """Indifference inside a 1-sigma box, euclidean distance outside."""

    ind = [o//2 for o in ord]
    WD = []
    for i in ind:
        if Js[i] > obs.Js[i] + obs.dJs[i]:
            WD.append((Js[i] - (obs.Js[i] + obs.dJs[i]))/obs.dJs[i])
        elif Js[i] < obs.Js[i] - obs.dJs[i]:
            WD.append((Js[i] - (obs.Js[i] - obs.dJs[i]))/obs.dJs[i])
        else:
            WD.append(0.0)
    WD = np.array(WD)
    return np.sqrt(sum(WD**2))

def mass(prof, obs):
    m = ah.mass_integral(*prof)
    return np.sqrt(((m - obs.M)/obs.dM)**2)

def k2(prof, obs):
    k2 = ah.lovek2(*prof)
    return np.sqrt(((k2 - obs.k2)/obs.dk2)**2)

def NMoI(I, obs):
    return np.sqrt(((I - obs.NMoI)/obs.dNMoI)**2)

# def period(P, obs):
#     return np.sqrt(((P - obs.P)/obs.dP)**2)

def rho0(prof, obs):
    """A prior on 1-bar density."""
    return np.sqrt(((prof[1][0] - obs.rho0)/(obs.drho0))**2)

def rhomax(prof, obs):
    """A penalty on unreasonably high central density."""

    romax = prof[1][-1]
    return (max(romax, obs.rhomax) - obs.rhomax)/obs.rhobar

### undocumented ad-hoc losses
def _french23(planet, obs):
    return np.sqrt((planet.s0 - obs.a0)**2)
