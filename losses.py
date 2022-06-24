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

def mass(tof, obs):
    m = ah.mass_integral(*tof)
    return np.sqrt(((m - obs.M)/obs.dM)**2)

def NMoI(I, obs):
    return np.sqrt(((I - obs.NMoI)/obs.dNMoI)**2)

def rhomax(tof, obs):
    """A penalty on unreasonably high central density."""

    romax = tof[1][-1]
    return (max(romax, obs.rhomax) - obs.rhomax)/obs.rhobar

def _mloss(prof,Js,obs,jflag):
    d2 = 0
    d2+= losses.euclid_Jnm(Js,obs,jflag)**2
    d2+= losses.mass(prof,obs)**2
    return np.sqrt(d2)
