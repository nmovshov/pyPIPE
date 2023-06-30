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
    return np.sqrt(((prof[1][0] - obs.rho0)/(obs.drho0))**2)

def rhomax(prof, obs):
    romax = prof[1][-1]
    return (max(romax, obs.rhomax) - obs.rhomax)/obs.rhobar

### undocumented ad-hoc losses
def _french23(planet, obs):
    def J2corr(J6p):
        J2= 3510.653e-6 + 0.40036 * (J6p-0.5e-6)
        return J2

    # Eq. C5, with \Delta a=0
    def J4corr(J6p):
        J4= -34.054e-6 + 1.0788 * (J6p-0.5e-6)
        return J4

    J2p, J4p, J6p = planet.Js[1:4]
    # terms in covariance matrix
    a = 0.389e-6**2
    b = 0.9859 * 0.389e-6 * 0.439e-6
    c = 0.439e-6**2
    # Eqs. C11, C12  eigenvalues
    lambda_1 = (a+c)/2 + np.sqrt(((a-c)/2)**2 + b**2)
    lambda_2 = (a+c)/2 - np.sqrt(((a-c)/2)**2 + b**2)
    # Eq. C17 - rotation of error ellipse in J2,J4 plane
    theta =np.arctan2(lambda_1 - a,b)
    # Eqs. C15, C16
    x = J2p - J2corr(J6p)
    y = J4p - J4corr(J6p)
    # Eqs. C13, C14
    xp =  x * np.cos(theta) + y * np.sin(theta)
    yp = -x * np.sin(theta) + y * np.cos(theta)
    # Eq. C6 - standard deviation of input J2', J4', J6' from occultation gravity solution
    r = np.sqrt(xp**2/lambda_1 + yp**2/lambda_2)
    return r
