#------------------------------------------------------------------------------
# Interior model of rotating fluid planet.
#------------------------------------------------------------------------------

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
        self.aos    = 1.0 # calculated equatorial to mean radius ratio (from tof<n>)
        self.G = 6.67430e-11; # m^3 kg^-1 s^-2 (2018 NIST reference)

        if obs is not None:
            self.set_observables(obs)

    def set_observables(self,obs):
        """ Copy physical properties from an observables struct."""
        self.mass = obs.M
        self.radius = obs.a0
        self.period = obs.P
        self.P0 = obs.P0