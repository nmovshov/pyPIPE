#------------------------------------------------------------------------------
# Interior model of rotating fluid planet.
#------------------------------------------------------------------------------
import numpy as np
import tof4
import tof7
from timeit import default_timer as timer

class TOFPlanet:
    """Interior model of rotating fluid planet.

    This class implements a model of a rotating fluid planet using Theory of
    Figures to calculate the hydrostatic equilibrium shape and resulting
    gravity field. A TOFPlanet object is defined by a densiy profile rho(s),
    supplied by the user and stored in the column vectors obj.si and
    obj.rhoi, indexed from the surface in. To complete the definition the
    user must also specify a mass, equatorial radius, and rotation period.
    With these a gravity field and equilibrium shape can be determined,
    with a call to obj.relax_to_HE(). Note, however, that the oblate shape
    calculated with relax_to_HE() preserves the mass and mean radius of the
    planet, but not the equatorial radius. A call to fix_radius()
    renormalizs the si vector to match the reference equatorial radius, at
    the cost of modifying the implied mass. A call to renormalize()
    modifies both si and rhoi to preserve the reference mass and equatorial
    radius, at the cost of modifying the assigned density. It is not
    possible to define mass, radius, and density simultaneously.

    """
    def __init__(self, obs=None):
        self.G = 6.67430e-11; # m^3 kg^-1 s^-2 (2018 NIST reference)
        self.opts = _default_opts() # holds user configurable options
        if obs is not None:
            self.set_observables(obs)

    def set_observables(self,obs):
        """ Copy physical properties from an observables struct."""
        self.mass = obs.M
        self.radius = obs.a0
        self.s0 = obs.s0
        self.P0 = obs.P0
        self.period = obs.P
        self.GM = self.G*self.mass
        self.wrot = 2*np.pi/self.period
        self.qrot = self.wrot**2*self.radius**3/self.GM
        self.mrot = self.wrot**2*self.s0**3/self.GM
        self.rhobar = self.mass/(4*np.pi/3*self.s0**3)

    def relax_to_HE(self, fixradius=True, fixmass=False, moi=False):
        """ Call tof<n> to obtain equilibrium shape and gravity."""

        if (self.opts['verbosity'] > 1):
            print('  Relaxing to hydrostatic equilibrium...')
        tic = timer()
        if self.opts['toforder'] == 4:
            tofun = tof4.tof4
        elif self.opts['toforder'] == 7:
            tofun = tof7.tof7
        else:
            raise ValueError('Unimplemented tof order')

        self.Js, out = tofun(self.si, self.rhoi, self.mrot,
            tol=self.opts['dJtol'], maxiter=self.opts['MaxIterHE'],
            xlevels=self.opts['xlevels'], calc_moi=moi)
        toc = timer() - tic

        if (self.opts['verbosity'] > 1):
            print('  Relaxing to hydrostatic equilibrium...done.')
            print(f' Elapsed time {toc:g} sec.')

        self.NMoI = out.NMoI
        self.ss = out.ss
        self.SS = out.SS
        self.A0 = out.A0
        self.aos = out.a0

        if fixradius:
            self.si = self.si*self.radius/(self.si[0]*self.aos)
            self.s0 = self.si[0]
            self.a0 = self.s0*self.aos
            self.M = _mass_int(self.si, self.rhoi)
        if fixmass:
            self.rhoi = self.rhoi*self.mass/self.M
            self.M = _mass_int(self.si, self.rhoi)
        if moi:
            self.NMoI = out.NMoI

def _mass_int(svec, dvec):
    """Trapz-integrate mass from rho(r) data."""
    from scipy.integrate import trapz
    return -4*np.pi*trapz(dvec*svec**2, x=svec)

def _default_opts():
    """Return options dict used by TOFPlanet class methods."""
    opts = {'toforder':4,
            'dJtol':1e-6,
            'MaxIterHE':60,
            'xlevels':-1,
            'verbosity':1
            }
    return opts
