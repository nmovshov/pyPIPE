"""Observed planetary values in consistent format.

The purpose of the observables module is to create a consistent and predictable
format for a structure containing observed values of a planet's vitals. The
typical usage is:

    from observables import <planet_name>

for the default, best-available data, or:

    from observables import <planet_name>_<mod_source>

for values and/or uncertainties modified to suit some purpose. For example,
observables.Saturn_tof4() modifies (increases) the gravity uncertainties to
match the estimated truncation error of 4th-order ToF.

The returned struct has the following fields:

  obs.pname            -  planet name
  obs.M, obs.dM        -  planet mass in kg, with uncertainty
  obs.a0               -  planet equatorial radius in meters
  obs.s0               -  planet mean radius in meters
  obs.P0               -  reference surface pressure on obs.a0, in Pa
  obs.T0, dT0          -  reference surface temperature, in K
  obs.rho0, drho0      -  reference density at (P0,T0), in kg/m^3
  obs.P, obs.dP        -  rotation period, in seconds
  obs.q, obs.dq        -  dimensionless rotation parameter, q=w^2*a0^3/(GM)
  obs.m, obs.dm        -  dimensionless rotation parameter, m=w^2*s0^3/(GM)
  obs.J<n>, obs.dJ<n>  -  n-th gravity harmonic and associated uncertainty.
                          J2-J14 fields are guaranteed to exists, although they
                          may contain NaN or Inf. The uncertainty value here
                          represents a convenient default 1-sigma value that we
                          commonly use, but is often adjusted on-the-fly in
                          client scripts. It sometimes is and sometimes isn't
                          the "official" uncertainty from the source dataset.
  obs.Js, obs.dJs      -  a vector of gravity harmonics and a vector of
                          corresponding default uncertainties. These are the
                          same values as in the individual obs.J<n> and
                          obs.dJ<n> fields; it's just sometimes more convenient
                          to use one or the other form. The vector forms are
                          always length 8, starting with J0 (always = -1) for
                          various reasons.

Important note about uncertainties: the d<x> quantities defined in the module
use reference values whose exact meaning may vary and may depend on context. It
is the user's job to decide if that value should be a 1-sigma, 2-sigma, or
uniform error bars, for example.
"""

import numpy as np

G = 6.67430e-11         # http://physics.nist.gov/cuu/index.html

class Jupiter:
    pname = 'jupiter'

    # Mass and radius, https://ssd.jpl.nasa.gov/ (2018)
    M  = 1898.187e24
    dM = 0.088e24
    a0 = 71492e3
    s0 = 69911e3

    # Boundary conditions
    P0 = 1e5            # The reference radius is the 1 bar level
    T0 = 165; dT0 = 5   # Lindal, G.F., 1992. Astrophys. J. 103, 967–982.
    rho0 = 0.169        # Protosolar ideal gas (mmw=2.319 amu) at (P0,T0)
    drho0 = 0.0051      # half the range of T0+/-dT0
    rhomax = 30000      # A guess, ANEOS serpentine at 50 Mbar is ~15000

    # Nominal rotation rate, https://ssd.jpl.nasa.gov/ (2018)
    P = 0.41354*24*3600
    w = 2*np.pi/P
    GM = G*M
    q = w**2*a0**3/GM
    m = w**2*s0**3/GM

    # Sometimes we use an estimate of roation period uncertainty
    dP = 30e-3 # Conservative! See e.g. Higgins et al. 1996
    dw = 2*np.pi/P**2*dP
    dq = 2*w*a0**3/GM*dw
    dm = 2*w*s0**3/GM*dw

    # Some methods (e.g. priors.py) use these fiducial values for scale...
    rhobar = M/(4*np.pi/3*s0**3)

    ## Gravity
    # Nominal coefficients from Iess et al. (2018) Table 1
    cfac = 71492e3/a0 # Iess et al. ref. radius converted to our equatorial radius
    J2  = 14696.572e-6*cfac**2
    J4  =  -586.609e-6*cfac**4
    J6  =    34.198e-6*cfac**6
    J8  =    -2.426e-6*cfac**8
    J10 =     0.172e-6*cfac**10
    J12 =     0.000e-6*cfac**12
    J14 =     0.000e-6*cfac**14

    # Formal uncertainties from Juno (we don't often use those)
    dJ2  = 0.014e-6*cfac**2
    dJ4  = 0.004e-6*cfac**4
    dJ6  = 0.009e-6*cfac**6
    dJ8  = 0.025e-6*cfac**8
    dJ10 = 0.069e-6*cfac**10
    dJ12 = np.inf
    dJ14 = np.inf

    # It is occasionally convenient to collect the Js and dJs in vectors.
    Js = np.array((-1, J2, J4, J6, J8, J10, J12, J14))
    dJs = np.array((0, dJ2, dJ4, dJ6, dJ8, dJ10, dJ12, dJ14))

    # A moment of inertia nominal value (not a real observation)
    NMoI = 0.2635
    dNMoI = 0.0005

    def __init__(self, **kwargs):
        """Customize instance."""

        # First override any defaults with directly supplied value
        for kw, val in kwargs.items():
            if kw in dir(self) and kwargs[kw] is not None:
                setattr(self, kw, val)

        # Let user override dJs and dM with more convenient *relative* sigs
        if 'J2_sig' in kwargs and kwargs['J2_sig'] is not None:
            self.dJ2 = abs(kwargs['J2_sig']*self.J2)
        if 'J4_sig' in kwargs and kwargs['J4_sig'] is not None:
            self.dJ4 = abs(kwargs['J4_sig']*self.J4)
        if 'J6_sig' in kwargs and kwargs['J6_sig'] is not None:
            self.dJ6 = abs(kwargs['J6_sig']*self.J6)
        if 'J8_sig' in kwargs and kwargs['J8_sig'] is not None:
            self.dJ8 = abs(kwargs['J8_sig']*self.J8)
        if 'J10_sig' in kwargs and kwargs['J10_sig'] is not None:
            self.dJ10 = abs(kwargs['J10_sig']*self.J10)
        if 'M_sig' in kwargs and kwargs['M_sig'] is not None:
            self.dM = kwargs['M_sig']*self.M

        # We have to manually reset Js and dJs for this instance
        self.Js = np.array(
                (-1,self.J2,self.J4,self.J6,self.J8,self.J10,self.J12,self.J14))
        self.dJs = np.array(
                (0,self.dJ2,self.dJ4,self.dJ6,self.dJ8,self.dJ10,self.dJ12,self.dJ14))

class Jupiter_winds(Jupiter):
    """Gravity nominals and uncertainties as in Miguel et al. (2022)."""
    J2  = Jupiter.J2  + 1.039e-6
    J4  = Jupiter.J4  - 0.076e-6
    J6  = Jupiter.J6  + 0.016e-6
    J8  = Jupiter.J8  + 0.053e-6
    J10 = Jupiter.J10 - 0.080e-6
    J12 = Jupiter.J12
    J14 = Jupiter.J14

    dJ2   = Jupiter.dJ2  + 0.354e-6
    dJ4   = Jupiter.dJ4  + 0.083e-6
    dJ6   = Jupiter.dJ6  + 0.076e-6
    dJ8   = Jupiter.dJ8  + 0.062e-6
    dJ10  = Jupiter.dJ10 + 0.042e-6
    dJ12  = Jupiter.dJ12
    dJ14  = Jupiter.dJ14

    Js = np.array([-1, J2, J4, J6, J8, J10, J12, J14])
    dJs = np.array([0, dJ2, dJ4, dJ6, dJ8, dJ10, dJ12, dJ14])

class Jupiter_tof4(Jupiter):
    """Modify gravity uncertainties to tof4 truncation error."""
    dJ2  = 1e-4*np.abs(Jupiter.J2)
    dJ4  = 3e-3*np.abs(Jupiter.J4)
    dJ6  = 3e-2*np.abs(Jupiter.J6)
    dJ8  = 3e-1*np.abs(Jupiter.J8)
    dJ10 = np.inf
    dJ12 = np.inf
    dJ14 = np.inf
    dJs = np.array((0, dJ2, dJ4, dJ6, dJ8, dJ10, dJ12, dJ14))

class Saturn:
    pname = 'saturn'

    # Mass and radius, https://ssd.jpl.nasa.gov/ (2018)
    M  = 568.336e24
    dM = 0.026e24
    a0 = 60268e3
    s0 = 58232e3

    # Boundary conditions
    P0 = 1e5            # The reference radius is the 1 bar level
    T0 = 140; dT0 = 4   # Nettelmann et al. (2013), Icarus 225, 548-557.
    rho0 = 0.199        # Protosolar ideal gas (mmw=2.319 amu) at (P0,T0)
    drho0 = 0.0057      # half the range of T0+/-dT0
    rhomax = 20000      # A guess, ANEOS serpentine at 50 Mbar is ~15000

    # Nominal rotation rate, Mankovich (2019) rounded to the minute
    P = 38040 # 10 hours 34 minutes
    w = 2*np.pi/P
    GM = G*M
    q = w**2*a0**3/GM
    m = w**2*s0**3/GM

    # Sometimes we use an estimate of roation period uncertainty
    dP = 120 # A 2-sigma ~= 2-minute spread of modern estimates
    dw = 2*np.pi/P**2*dP
    dq = 2*w*a0**3/GM*dw
    dm = 2*w*s0**3/GM*dw

    # Some methods (e.g. priors.py) use these fiducial values for scale...
    rhobar = M/(4*np.pi/3*s0**3)

    ## Gravity
    # Nominal coefficients from Iess et al. (2019)
    cfac = 60330e3/a0 # Iess et al. ref. radius converted to our equatorial radius
    J2  = +16290.573e-6*cfac**2
    J4  =   -935.314e-6*cfac**4
    J6  =    +86.340e-6*cfac**6
    J8  =    -14.624e-6*cfac**8
    J10 =     +4.672e-6*cfac**10
    J12 =     -0.000e-6*cfac**12
    J14 =     +0.000e-6*cfac**14

    # Formal uncertainties from Cassini (we don't often use those)
    dJ2  = 0.028e-6*cfac**2
    dJ4  = 0.037e-6*cfac**4
    dJ6  = 0.087e-6*cfac**6
    dJ8  = 0.205e-6*cfac**8
    dJ10 = 0.420e-6*cfac**10
    dJ12 = np.inf
    dJ14 = np.inf

    # It is occasionally convenient to collect the Js and dJs in vectors.
    Js = np.array((-1, J2, J4, J6, J8, J10, J12, J14))
    dJs = np.array((0, dJ2, dJ4, dJ6, dJ8, dJ10, dJ12, dJ14))

    def __init__(self, **kwargs):
        """Customize instance."""

        # First override any defaults with directly supplied value
        for kw, val in kwargs.items():
            if kw in dir(self) and kwargs[kw] is not None:
                setattr(self, kw, val)

        # Let user override dJs and dM with more convenient *relative* sigs
        if 'J2_sig' in kwargs and kwargs['J2_sig'] is not None:
            self.dJ2 = abs(kwargs['J2_sig']*self.J2)
        if 'J4_sig' in kwargs and kwargs['J4_sig'] is not None:
            self.dJ4 = abs(kwargs['J4_sig']*self.J4)
        if 'J6_sig' in kwargs and kwargs['J6_sig'] is not None:
            self.dJ6 = abs(kwargs['J6_sig']*self.J6)
        if 'J8_sig' in kwargs and kwargs['J8_sig'] is not None:
            self.dJ8 = abs(kwargs['J8_sig']*self.J8)
        if 'J10_sig' in kwargs and kwargs['J10_sig'] is not None:
            self.dJ10 = abs(kwargs['J10_sig']*self.J10)
        if 'M_sig' in kwargs and kwargs['M_sig'] is not None:
            self.dM = kwargs['M_sig']*self.M

        # We have to manually reset Js and dJs for this instance
        self.Js = np.array(
                (-1,self.J2,self.J4,self.J6,self.J8,self.J10,self.J12,self.J14))
        self.dJs = np.array(
                (0,self.dJ2,self.dJ4,self.dJ6,self.dJ8,self.dJ10,self.dJ12,self.dJ14))

class Saturn_tof4(Saturn):
    """Modify gravity uncertainties to tof4 truncation error."""
    dJ2  = 1e-4*np.abs(Saturn.J2)
    dJ4  = 3e-3*np.abs(Saturn.J4)
    dJ6  = 3e-2*np.abs(Saturn.J6)
    dJ8  = 3e-1*np.abs(Saturn.J8)
    dJ10 = np.inf
    dJ12 = np.inf
    dJ14 = np.inf
    dJs = np.array((0, dJ2, dJ4, dJ6, dJ8, dJ10, dJ12, dJ14))

class Saturn_winds(Saturn):
    """Gravity uncertainties reflecting potential deep wind contribution.

    See fig. 4 in Galanti, E., & Kaspi, Y. (2017). The Astrophysical Journal,
    843(2), L25.
    """
    dJ2  = 15e-6 # I interpret Galanti+2017 (fig. 4) as 2-sigma=3e-5
    dJ4  = 5e-6  # Represents *generous* ToF model error + deep winds
    dJ6  = 5e-6  # Represents *generous* ToF model error + deep winds
    dJ8  = np.inf
    dJ10 = np.inf
    dJ12 = np.inf
    dJ14 = np.inf
    dJs = np.array((0, dJ2, dJ4, dJ6, dJ8, dJ10, dJ12, dJ14))

class Uranus:
    pname = 'uranus'

    # Mass and radius, https://ssd.jpl.nasa.gov/ (2018)
    M = 86.8127e24
    dM = 0.004e24
    a0 = 25559e3
    s0 = 25362e3

    # Boundary conditions
    P0 = 1e5            # The reference radius is the 1 bar level
    T0 = 76; dT0 = 2    # Lindal, G.F., 1992. Astrophys. J. 103, 967–982.
    rho0 = 0.367        # Protosolar ideal gas (mmw=2.319 amu) at (P0,T0)
    drho0 = 0.0097      # half the range of T0+/-dT0
    rhomax = 20000      # Generous guess

    # Nominal rotation rate, https://ssd.jpl.nasa.gov/ (2018)
    P = 0.71833*24*3600
    w = 2*np.pi/P
    GM = G*M
    q = w**2*a0**3/GM
    m = w**2*s0**3/GM

    # Sometimes we use an estimate of roation period uncertainty
    dP = 600 # Basically a wild guess, see e.g. Podolak and Helled 2012
    dw = 2*np.pi/P**2*dP
    dq = 2*w*a0**3/GM*dw
    dm = 2*w*s0**3/GM*dw

    # Some methods (e.g. priors.py) use these fiducial values for scale...
    rhobar = M/(4*np.pi/3*s0**3)

    ## Gravity
    # Nominal coefficients from Jacobson (2014) table 12
    cfac = 25559e3/a0 # Jacobson (2014) radius converted to our equatorial radius
    J2  = +3510.7e-6*cfac**2
    J4  =   -34.2e-6*cfac**4
    J6  =          0*cfac**6
    J8  =          0*cfac**8
    J10 =          0*cfac**10
    J12 =          0*cfac**12
    J14 =          0*cfac**14

    # Recommended uncertainties
    dJ2  = 0.7e-6*cfac**2
    dJ4  = 1.3e-6*cfac**4
    dJ6  = np.inf
    dJ8  = np.inf
    dJ10 = np.inf
    dJ12 = np.inf
    dJ14 = np.inf

    # It is occasionally convenient to collect the Js and dJs in vectors.
    Js = np.array((-1, J2, J4, J6, J8, J10, J12, J14))
    dJs = np.array((0, dJ2, dJ4, dJ6, dJ8, dJ10, dJ12, dJ14))

    def __init__(self, **kwargs):
        """Customize instance."""

        # First override any defaults with directly supplied value
        for kw, val in kwargs.items():
            if kw in dir(self) and kwargs[kw] is not None:
                setattr(self, kw, val)

        # Let user override dJs and dM with more convenient *relative* sigs
        if 'J2_sig' in kwargs and kwargs['J2_sig'] is not None:
            self.dJ2 = abs(kwargs['J2_sig']*self.J2)
        if 'J4_sig' in kwargs and kwargs['J4_sig'] is not None:
            self.dJ4 = abs(kwargs['J4_sig']*self.J4)
        if 'J6_sig' in kwargs and kwargs['J6_sig'] is not None:
            self.dJ6 = abs(kwargs['J6_sig']*self.J6)
        if 'J8_sig' in kwargs and kwargs['J8_sig'] is not None:
            self.dJ8 = abs(kwargs['J8_sig']*self.J8)
        if 'J10_sig' in kwargs and kwargs['J10_sig'] is not None:
            self.dJ10 = abs(kwargs['J10_sig']*self.J10)
        if 'M_sig' in kwargs and kwargs['M_sig'] is not None:
            self.dM = kwargs['M_sig']*self.M

        # We have to manually reset Js and dJs for this instance
        self.Js = np.array(
                (-1,self.J2,self.J4,self.J6,self.J8,self.J10,self.J12,self.J14))
        self.dJs = np.array(
                (0,self.dJ2,self.dJ4,self.dJ6,self.dJ8,self.dJ10,self.dJ12,self.dJ14))

class Uranus_ppwd(Uranus):
    J2  =  Uranus.J2
    J4  =  Uranus.J4
    J6  =  5.1769e-7
    J8  = -1.0421e-8
    J10 =  2.5672e-10
    J12 = -7.2879e-12
    J14 =  2.4279e-13

    dJ2  = 1e-6*J2
    dJ4  = 1e-5*J4
    dJ6  = 1e-4*J6
    dJ8  = 1e-4*J8
    dJ10 = 1e-2*J10
    dJ12 = 1e-0*J12
    dJ14 = 1e-0*J14

    Js = np.array((-1, J2, J4, J6, J8, J10, J12, J14))
    dJs = np.array((0, dJ2, dJ4, dJ6, dJ8, dJ10, dJ12, dJ14))

class Uranus_uncertain_rotation(Uranus_ppwd):
    dP = 1800
    dw = 2*np.pi/Uranus.P**2*dP
    dq = 2*Uranus.w*Uranus.a0**3/Uranus.GM*dw
    dm = 2*Uranus.w*Uranus.s0**3/Uranus.GM*dw

class Neptune:
    pname = 'neptune'

    # Mass and radius, https://ssd.jpl.nasa.gov/ (2018)
    M = 102.4126e24
    dM = 0.0048e24
    a0 = 24764e3
    s0 = 24622e3

    # Boundary conditions
    P0 = 1e5            # The reference radius is the 1 bar level
    T0 = 72; dT0 = 2    # Lindal, G.F., 1992. Astrophys. J. 103, 967–982.
    rho0 = 0.387        # Protosolar ideal gas (mmw=2.319 amu) at (P0,T0)
    drho0 = 0.0108      # half the range of T0+/-dT0
    rhomax = 20000      # Generous guess

    # Nominal rotation rate, https://ssd.jpl.nasa.gov/ (2018)
    P = 0.67125*24*3600
    w = 2*np.pi/P
    GM = G*M
    q = w**2*a0**3/GM
    m = w**2*s0**3/GM

    # Sometimes we use an estimate of roation period uncertainty
    dP = 600 # Basically a wild guess, see e.g. Podolak and Helled 2012
    dw = 2*np.pi/P**2*dP
    dq = 2*w*a0**3/GM*dw
    dm = 2*w*s0**3/GM*dw

    # Some methods (e.g. priors.py) use these fiducial values for scale...
    rhobar = M/(4*np.pi/3*s0**3)

    ## Gravity
    # Nominal coefficients from Jacobson (2009) table 5
    cfac = 25225e3/a0 # Jacobson (2009) radius converted to our equatorial radius
    J2  = +3408.4e-6*cfac**2
    J4  =   -33.4e-6*cfac**4
    J6  =          0*cfac**6
    J8  =          0*cfac**8
    J10 =          0*cfac**10
    J12 =          0*cfac**12
    J14 =          0*cfac**14

    # Recommended uncertainties
    dJ2  = 4.5e-6*cfac**2
    dJ4  = 2.9e-6*cfac**4
    dJ6  = np.inf
    dJ8  = np.inf
    dJ10 = np.inf
    dJ12 = np.inf
    dJ14 = np.inf

    # It is occasionally convenient to collect the Js and dJs in vectors.
    Js = np.array((-1, J2, J4, J6, J8, J10, J12, J14))
    dJs = np.array((0, dJ2, dJ4, dJ6, dJ8, dJ10, dJ12, dJ14))

    def __init__(self, **kwargs):
        """Customize instance."""

        # First override any defaults with directly supplied value
        for kw, val in kwargs.items():
            if kw in dir(self) and kwargs[kw] is not None:
                setattr(self, kw, val)

        # Let user override dJs and dM with more convenient *relative* sigs
        if 'J2_sig' in kwargs and kwargs['J2_sig'] is not None:
            self.dJ2 = abs(kwargs['J2_sig']*self.J2)
        if 'J4_sig' in kwargs and kwargs['J4_sig'] is not None:
            self.dJ4 = abs(kwargs['J4_sig']*self.J4)
        if 'J6_sig' in kwargs and kwargs['J6_sig'] is not None:
            self.dJ6 = abs(kwargs['J6_sig']*self.J6)
        if 'J8_sig' in kwargs and kwargs['J8_sig'] is not None:
            self.dJ8 = abs(kwargs['J8_sig']*self.J8)
        if 'J10_sig' in kwargs and kwargs['J10_sig'] is not None:
            self.dJ10 = abs(kwargs['J10_sig']*self.J10)
        if 'M_sig' in kwargs and kwargs['M_sig'] is not None:
            self.dM = kwargs['M_sig']*self.M

        # We have to manually reset Js and dJs for this instance
        self.Js = np.array(
                (-1,self.J2,self.J4,self.J6,self.J8,self.J10,self.J12,self.J14))
        self.dJs = np.array(
                (0,self.dJ2,self.dJ4,self.dJ6,self.dJ8,self.dJ10,self.dJ12,self.dJ14))

class Neptune_ppwd(Neptune):
    J2  =  Neptune.J2
    J4  =  Neptune.J4
    J6  =  5.7563e-7
    J8  = -1.2049e-8
    J10 =  2.9915e-10
    J12 = -8.2982e-12
    J14 =  2.6246e-13

    dJ2  = 1e-6*J2
    dJ4  = 1e-5*J4
    dJ6  = 1e-4*J6
    dJ8  = 1e-4*J8
    dJ10 = 1e-2*J10
    dJ12 = 1e-0*J12
    dJ14 = 1e-0*J14

    Js = np.array((-1, J2, J4, J6, J8, J10, J12, J14))
    dJs = np.array((0, dJ2, dJ4, dJ6, dJ8, dJ10, dJ12, dJ14))

class Neptune_uncertain_rotation(Neptune_ppwd):
    dP = 1800
    dw = 2*np.pi/Neptune.P**2*dP
    dq = 2*Neptune.w*Neptune.a0**3/Neptune.GM*dw
    dm = 2*Neptune.w*Neptune.s0**3/Neptune.GM*dw

class HAT_P_13_b:
    pname = 'HAT-P-13-b'

    # Bakos et al. 2009
    M  = 0.85*Jupiter().M
    dM = 0.05*Jupiter().M
    a0 = 1.28*Jupiter().a0
    s0 = a0

    # From Batygin et al. 2009
    Qlo = 10000
    Qhi = 300000
    k2 = 0.2705
    dk2 = 0.155

    # Boundary conditions (Bakos et al. 2009)
    P0 = 1e5            # The reference radius is the 1 bar level
    T0 = 1653; dT0 = 45 # Table 4
    rho0 = 0.0169       # Protosolar ideal gas (mmw=2.319 amu) at (P0,T0)
    drho0 = 0.0005      # half the range of T0+/-dT0
    rhomax = 30000      # A guess, ANEOS serpentine at 50 Mbar is ~15000

    # Nominal rotation rate, https://ssd.jpl.nasa.gov/ (2018)
    P = np.inf
    w = 2*np.pi/P
    GM = G*M
    q = w**2*a0**3/GM
    m = w**2*s0**3/GM

    # Some methods (e.g. priors.py) use these fiducial values for scale...
    rhobar = None

    # Empty placeholders for Js and dJs
    Js = np.array((-1., 0, 0, 0, 0, 0, 0, 0))
    dJs = np.array((0,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf))

    def __init__(self, **kwargs):
        """Customize instance."""

        # Manual override any defaults with directly supplied value
        for kw, val in kwargs.items():
            if kw in dir(self) and kwargs[kw] is not None:
                setattr(self, kw, val)

if __name__ == "__main__":
    print(Saturn)
