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
        self.mass   = 0.0 # reference mass
        self.radius = 0.0 # reference radius (equatorial!)
        self.period = 0.0 # reference rotation period
        self.P0     = 0.0 # reference pressure
        self.si     = 0.0 # vector of mean radii (top down, s0=si[0] is outer radius)
        self.rhoi   = 0.0 # vector of densities on si grid
        self.Js     = 0.0 # external gravity coefficients (returned by tof<n>)
        self.M      = 0.0 # calculated mass
        self.a0     = 0.0 # calculated equatorial radius
        self.s0     = 0.0 # surface mean radius (another name for si[0]
        self.mi     = 0.0 # cumulative mass below si
        self.ai     = 0.0 # equatorial radii on level surfaces
        self.rhobar = 0.0 # calculated mean density
        self.wrot   = 0.0 # rotation frequency, 2pi/period
        self.qrot   = 0.0 # rotation parameter wrot^2a0^3/GM
        self.mrot   = 0.0 # rotation parameter, wrot^2s0^3/GM
        self.ss     = 0.0 # shape functions (returned by tof<n>)
        self.SS     = 0.0 # shape functions (returned by tof<n>)
        self.A0     = 0.0 # dimensionless potential (returned by tof<n>)
        self.Js     = 0.0 # external gravity coefficients (returned by tof<n>)
        self.GM     = 0.0 # mass parameter
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

    def relax_to_HE(self):
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
            xlevels=self.opts['xlevels'])
        toc = timer() - tic

        self.aos = out.a0
        # self.ss = structfun(@flipud, out.ss, 'UniformOutput', false);
        # self.SS = structfun(@flipud, out.SS, 'UniformOutput', false);
        # self.A0 = flipud(out.A0)

        if (self.opts['verbosity'] > 1):
            print('  Relaxing to hydrostatic equilibrium...done.')
            print(f' Elapsed time {toc:g} sec.')

def _default_opts():
    """Return options dict used by TOFPlanet class methods."""
    opts = {'toforder':4,
            'dJtol':1e-6,
            'drhotol':1e-6,
            'MaxIterBar':60,
            'MaxIterHE':60,
            'xlevels':-1,
            'verbosity':1
            }
    return opts
