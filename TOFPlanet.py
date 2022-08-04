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
    def __init__(self, obs=None, name='planet'):
        self.G = 6.67430e-11; # m^3 kg^-1 s^-2 (2018 NIST reference)
        self.opts = _default_opts() # holds user configurable options
        if obs is None:
            obs = _default_planet()
        self.set_observables(obs)
        self.name = name

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

    def relax_to_HE(self, fixradius=True, fixmass=False,
                    moi=False, pressure=False):
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
        self.A0 = np.flip(out.A0)
        self.aos = out.a0

        if fixradius:
            self.si = self.si*self.radius/(self.si[0]*self.aos)
            self.s0 = self.si[0]
            self.a0 = self.s0*self.aos
            self.M = _mass_int(self.si, self.rhoi)
        if fixmass:
            self.rhoi = self.rhoi*self.mass/self.M
            self.M = _mass_int(self.si, self.rhoi)

        if pressure:
            r = self.si
            rho = self.rhoi
            U = -self.G*self.mass/self.s0**3*self.si**2*self.A0
            gradU = np.zeros_like(r)
            gradU[0] = (U[0] - U[1])/(r[0] - r[1])
            gradU[1:-1] = (U[0:-2] - U[2:])/(r[0:-2] - r[2:])
            gradU[-1] = (U[-2] - U[-1])/(r[-2] - r[-1])
            intgrnd = rho*gradU
            P = np.zeros_like(r)
            P[0] = self.P0
            for k in range(P.size-1): # note integrate downward
                P[k+1] = P[k] + 0.5*(r[k] - r[k+1])*(intgrnd[k] + intgrnd[k+1])
            self.Pi = P

    ### Visualizers
    def plot_rho_of_r(self):
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(8,6))
        x = np.append(self.si/self.s0, 0)
        y = np.append(self.rhoi, self.rhoi[-1])/1000
        plt.plot(x, y, lw=2, label=self.name)
        plt.xlabel(r'Level surface radius, $s/s_0$', fontsize=12)
        plt.ylabel(r'$\rho$ [1000 kg/m$^3$]', fontsize=12)
        plt.show(block=False)

    def plot_P_of_r(self):
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(8,6))
        x = self.si/self.s0
        y = self.Pi/1e9
        plt.plot(x, y, lw=2, label=self.name)
        plt.xlabel(r'Level surface radius, $s/s_0$', fontsize=12)
        plt.ylabel(r'$P$ [GPa]', fontsize=12)
        plt.show(block=False)

    def plot_P_of_rho(self):
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(8,6))
        x = self.rhoi/1000
        y = self.Pi/1e9
        plt.plot(x, y, lw=2, label=self.name)
        plt.xlabel(r'$\rho$ [1000 kg/m$^3$]', fontsize=12)
        plt.ylabel(r'$P$ [GPa]', fontsize=12)
        plt.show(block=False)

    def plot_rho_of_P(self):
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(8,6))
        x = self.Pi/1e9
        y = self.rhoi/1000
        plt.loglog(x, y, lw=2, label=self.name)
        plt.xlabel(r'$P$ [GPa]', fontsize=12)
        plt.ylabel(r'$\rho$ [1000 kg/m$^3$]', fontsize=12)
        plt.xlim(left=1e-3)
        plt.show(block=False)

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

class _default_planet:
    """Use this to prefill critical TP fields with reasonable values."""
    M  = 1898.187e24
    a0 = 71492e3
    s0 = 69911e3
    P0 = 1e5
    P = 0.41354*24*3600

def _test():
    obs = _default_planet()
    tp = TOFPlanet(obs)
    N = 4096
    tp.si = np.linspace(1, 1/N, N)
    a = - 15*obs.M/8/np.pi/obs.s0**3
    tp.rhoi = a*tp.si**2 - a
    tp.relax_to_HE(moi=True, pressure=True)
    print("J0 = {}".format(tp.Js[0]))
    print("J2 = {}".format(tp.Js[1]))
    print("J4 = {}".format(tp.Js[2]))
    print("J6 = {}".format(tp.Js[3]))
    print("J8 = {}".format(tp.Js[4]))
    print("I = {}".format(tp.NMoI))
    print("P_center = {}".format(tp.Pi[-1]))
    print("")

if __name__ == '__main__':
    _test()
