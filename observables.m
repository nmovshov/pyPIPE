classdef observables
%OBSERVABLES Observed planetary values in consistent format.
% The purpose of the observables class is to create a consistent and
% predictable format for a structure containing observed values of a planet's
% vitals. The typical usage is:
%
%   obs = observables.<planet_name>()
%
% for the default, best-available data, or:
%
%   obs = observables.<planet_name>_<mod_source>()
%
% for values and/or uncertainties modified to suit some purpose. For example,
% observables.Saturn4() modifies (increases) the gravity uncertainties to
% match the estimated truncation error of 4th-order ToF.
%
% The returned struct has the following fields:
%
%   obs.pname            -  planet name
%   obs.M, obs.dM        -  planet mass in kg, with uncertainty
%   obs.a0               -  planet equatorial radius in meters
%   obs.s0               -  planet mean radius in meters
%   obs.P0               -  reference surface pressure on obs.a0, in Pa
%   obs.T0, dT0          -  reference surface temperature, in K
%   obs.rho0, drho0      -  reference density at (P0,T0), in kg/m^3
%   obs.P, obs.dP        -  rotation period, in seconds
%   obs.q, obs.dq        -  dimensionless rotation parameter, q=w^2*a0^3/(GM)
%   obs.m, obs.dm        -  dimensionless rotation parameter, m=w^2*s0^3/(GM)
%   obs.J<n>, obs.dJ<n>  -  n-th gravity harmonic and associated uncertainty.
%                           J2-J14 fields are guaranteed to exists, although
%                           they may contain NaN or Inf. The uncertainty value
%                           here represents a convenient default 1-sigma value
%                           that we commonly use, but is often adjusted
%                           on-the-fly in client scripts. It sometimes is and
%                           sometimes isn't the "official" uncertainty from the
%                           source dataset.
%   obs.Js, obs.dJs      -  a vector of gravity harmonics and a vector of
%                           corresponding default uncertainties. These are the
%                           same values as in the individual obs.J<n> and
%                           obs.dJ<n> fields; it's just sometimes more
%                           convenient to use one or the other form. The vector
%                           forms are always length 8, starting with J0 (always
%                           = -1) for various reasons.
%
% Important note about uncertainties: the d<x> quantities defined in the module
% use reference values whose exact meaning may vary and may depend on context.
% It is the user's job to decide if that value should be a 1-sigma, 2-sigma, or
% uniform error bars, for example.

properties
    pname
    M
    dM
    a0
    s0
    rhobar
    P0
    T0
    dT0
    rho0
    drho0
    rhomax
    P
    GM
    w
    q
    m
    dP
    dw
    dq
    dm
    J2
    J4
    J6
    J8
    J10
    J12
    J14
    dJ2
    dJ4
    dJ6
    dJ8
    dJ10
    dJ12
    dJ14
    Js
    dJs
    NMoI
    dNMoI
end
properties (GetAccess=private)
    G = 6.67430e-11 % http://physics.nist.gov/cuu/index.html
end

%% Constructor
methods
    function obj = observables(varargin)
        %OBSERVABLES Construct an instance of the base class.
        assert(mod(nargin,2) == 0, "Odd number of name/value pairs.")
        while ~isempty(varargin)
            name = varargin{end-1};
            value = varargin{end};
            varargin(end-1:end) = [];
            if isprop(obj,name)
                obj.(name) = value;
            else
                warning("Unknown property: %s", name)
            end
        end
    end
end % constructor

%% Planets (as static methods)
methods (Static)
    function obs = Jupiter()
        obs = observables();
        obs.pname = 'jupiter';

        % Mass and radius, https://ssd.jpl.nasa.gov/ (2018)
        obs.M  = 1898.187e24;
        obs.dM = 0.088e24;
        obs.a0 = 71492e3;
        obs.s0 = 69911e3;

        % Boundary conditions
        obs.P0 = 1e5;              % The reference radius is the 1 bar level
        obs.T0 = 165; obs.dT0 = 5; % Lindal, G.F., 1992. Astrophys. J. 103, 967–982.
        obs.rho0 = 0.169;          % Protosolar ideal gas (mmw=2.319 amu) at (P0,T0)
        obs.drho0 = 0.0051;        % half the range of T0+/-dT0
        obs.rhomax = 30000;        % A guess, ANEOS serpentine at 50 Mbar is ~15000

        % Nominal rotation rate, https://ssd.jpl.nasa.gov/ (2018)
        obs.P = 0.41354*24*3600;
        obs.w = 2*pi/obs.P;
        obs.GM = obs.G*obs.M;
        obs.q = obs.w^2*obs.a0^3/obs.GM;
        obs.m = obs.w^2*obs.s0^3/obs.GM;

        % Sometimes we use an estimate of roation period uncertainty
        obs.dP = 30e-3; % Conservative! See e.g. Higgins et al. 1996
        obs.dw = 2*pi/obs.P^2*obs.dP;
        obs.dq = 2*obs.w*obs.a0^3/obs.GM*obs.dw;
        obs.dm = 2*obs.w*obs.s0^3/obs.GM*obs.dw;

        % Some methods (e.g. priors.py) use these fiducial values for scale...
        obs.rhobar = obs.M/(4*pi/3*obs.s0^3);

        % Gravity
        % Nominal coefficients from Iess et al. (2018) Table 1
        cfac = 71492e3/obs.a0; % Iess et al. ref. radius to our equatorial radius
        obs.J2  = 14696.572e-6*cfac^2;
        obs.J4  =  -586.609e-6*cfac^4;
        obs.J6  =    34.198e-6*cfac^6;
        obs.J8  =    -2.426e-6*cfac^8;
        obs.J10 =     0.172e-6*cfac^10;
        obs.J12 =     0.000e-6*cfac^12;
        obs.J14 =     0.000e-6*cfac^14;
        % Formal uncertainties from Juno (we don't often use those)
        obs.dJ2  = 0.014e-6*cfac^2;
        obs.dJ4  = 0.004e-6*cfac^4;
        obs.dJ6  = 0.009e-6*cfac^6;
        obs.dJ8  = 0.025e-6*cfac^8;
        obs.dJ10 = 0.069e-6*cfac^10;
        obs.dJ12 = inf;
        obs.dJ14 = inf;
        % It is occasionally convenient to collect the Js and dJs in vectors.
        obs.Js = [-1, obs.J2, obs.J4, obs.J6, obs.J8, obs.J10, obs.J12, obs.J14];
        obs.dJs = [0, obs.dJ2, obs.dJ4, obs.dJ6, obs.dJ8, obs.dJ10, obs.dJ12, obs.dJ14];

        % A moment of inertia nominal value (not a real observation)
        obs.NMoI = 0.2635;
        obs.dNMoI = 0.0005;
    end

    function obs = Jupiter_tof4()
        % Modify gravity uncertainties to tof4 truncation error.
        obs = observables.Jupiter();
        obs.dJ2  = 1e-4*abs(obs.J2);
        obs.dJ4  = 3e-3*abs(obs.J4);
        obs.dJ6  = 3e-2*abs(obs.J6);
        obs.dJ8  = 3e-1*abs(obs.J8);
        obs.dJ10 = inf;
        obs.dJ12 = inf;
        obs.dJ14 = inf;
        obs.dJs = [0, obs.dJ2, obs.dJ4, obs.dJ6, obs.dJ8,...
                      obs.dJ10, obs.dJ12, obs.dJ14];
    end

    function obs = Jupiter_winds()
        % Gravity nominals and uncertainties as in Miguel et al. (2022)
        obs = observables.Jupiter();
        obs.J2  = obs.J2  + 1.039e-6;
        obs.J4  = obs.J4  - 0.076e-6;
        obs.J6  = obs.J6  + 0.016e-6;
        obs.J8  = obs.J8  + 0.053e-6;
        obs.J10 = obs.J10 - 0.080e-6;

        obs.dJ2   = obs.dJ2  + 0.354e-6;
        obs.dJ4   = obs.dJ4  + 0.083e-6;
        obs.dJ6   = obs.dJ6  + 0.076e-6;
        obs.dJ8   = obs.dJ8  + 0.062e-6;
        obs.dJ10  = obs.dJ10 + 0.042e-6;
        
        obs.Js = [-1, obs.J2, obs.J4, obs.J6, obs.J8, obs.J10, obs.J12, obs.J14];
        obs.dJs = [0, obs.dJ2, obs.dJ4, obs.dJ6, obs.dJ8,...
                      obs.dJ10, obs.dJ12, obs.dJ14];
    end

    function obs = Saturn()
        obs = observables();
        obs.pname = 'Saturn';

        % Mass and radius, https://ssd.jpl.nasa.gov/ (2018)
        obs.M  = 568.336e24;
        obs.dM = 0.026e24;
        obs.a0 = 60268e3;
        obs.s0 = 58232e3;

        % Boundary conditions
        obs.P0 = 1e5;              % The reference radius is the 1 bar level
        obs.T0 = 140; obs.dT0 = 4; % Nettelmann et al. (2013)
        obs.rho0 = 0.199;          % Ideal gas (mmw=2.319 amu) at (P0,T0)
        obs.drho0 = 0.0057;        % half the range of T0+/-dT0
        obs.rhomax = 20000;        % ANEOS serpentine at 50 Mbar is ~15000

        % Nominal rotation rate e.g. Mankovich (2019)
        obs.P = 38014 % 10 hours 33 minutes 34 seconds
        obs.w = 2*pi/obs.P;
        obs.GM = obs.G*obs.M;
        obs.q = obs.w^2*obs.a0^3/obs.GM;
        obs.m = obs.w^2*obs.s0^3/obs.GM;

        % Sometimes we use an estimate of roation period uncertainty
        obs.dP = 55; % A 2-sigma ~= 1-minute spread of modern estimates
        obs.dw = 2*pi/obs.P^2*obs.dP;
        obs.dq = 2*obs.w*obs.a0^3/obs.GM*obs.dw;
        obs.dm = 2*obs.w*obs.s0^3/obs.GM*obs.dw;

        % Some methods (e.g. priors.py) use these fiducial values for scale...
        obs.rhobar = obs.M/(4*pi/3*obs.s0^3);

        % Gravity
        % Nominal coefficients from Iess et al. (2019)
        cfac = 60330e3/obs.a0; % Iess et al. ref. radius to equatorial radius
        obs.J2  = +16290.573e-6*cfac^2;
        obs.J4  =   -935.314e-6*cfac^4;
        obs.J6  =    +86.340e-6*cfac^6;
        obs.J8  =    -14.624e-6*cfac^8;
        obs.J10 =     +4.672e-6*cfac^10;
        obs.J12 =     -0.000e-6*cfac^12;
        obs.J14 =     +0.000e-6*cfac^14;

        % Formal uncertainties from Cassini (we don't often use those)
        obs.dJ2  = 0.028e-6*cfac^2;
        obs.dJ4  = 0.037e-6*cfac^4;
        obs.dJ6  = 0.087e-6*cfac^6;
        obs.dJ8  = 0.205e-6*cfac^8;
        obs.dJ10 = 0.420e-6*cfac^10;
        obs.dJ12 = inf;
        obs.dJ14 = inf;

        % It is occasionally convenient to collect the Js and dJs in vectors.
        obs.Js = [-1,obs.J2,obs.J4,obs.J6,obs.J8,obs.J10,obs.J12,obs.J14];
        obs.dJs = [0,obs.dJ2,obs.dJ4,obs.dJ6,obs.dJ8,...
                     obs.dJ10,obs.dJ12,obs.dJ14];
    end

    function obs = Saturn1()
        % Slightly relaxed gravity uncertainties, rounded up in oom..
        obs = observables.Saturn();
        obs.dJ2  = 1e-7;
        obs.dJ4  = 1e-7;
        obs.dJ6  = 1e-7;
        obs.dJ8  = inf;
        obs.dJ10 = inf;
        obs.dJ12 = inf;
        obs.dJ14 = inf;
        obs.dJs = [0, obs.dJ2, obs.dJ4, obs.dJ6, obs.dJ8,...
                      obs.dJ10, obs.dJ12, obs.dJ14];
    end

    function obs = Saturn4()
        % Modify gravity uncertainties to tof4 truncation error.
        obs = observables.Saturn();
        obs.dJ2  = 1e-4*abs(obs.J2);
        obs.dJ4  = 3e-3*abs(obs.J4);
        obs.dJ6  = 3e-2*abs(obs.J6);
        obs.dJ8  = 3e-1*abs(obs.J8);
        obs.dJ10 = inf;
        obs.dJ12 = inf;
        obs.dJ14 = inf;
        obs.dJs = [0, obs.dJ2, obs.dJ4, obs.dJ6, obs.dJ8,...
                      obs.dJ10, obs.dJ12, obs.dJ14];
    end

    function obs = Saturn9()
        % Gravity uncertainties reflecting potential deep wind contribution.
        %
        % See fig. 4 in Galanti, E., & Kaspi, Y. (2017). The Astrophysical
        % Journal, 843(2), L25.
        obs = observables.Saturn();
        obs.dJ2  = 15e-6; % I interpret Galanti+2017 (fig. 4) as 2-sigma=3e-5
        obs.dJ4  = 5e-6;  % Represents *generous* ToF model error + deep winds
        obs.dJ6  = 5e-6;  % Represents *generous* ToF model error + deep winds
        obs.dJ8  = inf;
        obs.dJ10 = inf;
        obs.dJ12 = inf;
        obs.dJ14 = inf;
        obs.dJs = [0,obs.dJ2,obs.dJ4,obs.dJ6,obs.dJ8,...
                     obs.dJ10,obs.dJ12,obs.dJ14];
    end

    function obs = Uranus()
        obs = observables();
        obs.pname = 'Uranus';

        % Mass and radius, https://ssd.jpl.nasa.gov/ (2018)
        obs.M = 86.8127e24;
        obs.dM = 0.004e24;
        obs.a0 = 25559e3;
        obs.s0 = 25362e3;

        % Boundary conditions
        obs.P0 = 1e5;             % The reference radius is the 1 bar level
        obs.T0 = 76; obs.dT0 = 2; % Lindal, G.F., 1992. Astrophys. J. 103
        obs.rho0 = 0.367;         % Ideal gas (mmw=2.319 amu) at (P0,T0)
        obs.drho0 = 0.0097;       % Half the range of T0+/-dT0
        obs.rhomax = 20000;       % Generous guess

        % Rotation rate, https://ssd.jpl.nasa.gov/ (2018)
        obs.P = 0.71833*24*3600;
        obs.w = 2*pi/obs.P;
        obs.GM = obs.G*obs.M;
        obs.q = obs.w^2*obs.a0^3/obs.GM;
        obs.m = obs.w^2*obs.s0^3/obs.GM;
        obs.dP = 600; % Basically a wild guess, e.g. Podolak and Helled 2012
        obs.dw = 2*pi/obs.P^2*obs.dP;
        obs.dq = 2*obs.w*obs.a0^3/obs.GM*obs.dw;
        obs.dm = 2*obs.w*obs.s0^3/obs.GM*obs.dw;

        % Some methods (e.g. priors.py) use these fiducial values for scale...
        obs.rhobar = obs.M/(4*pi/3*obs.s0^3);

        % Gravity
        % Nominal coefficients from Jacobson (2014) table 12
        cfac = 25559e3/obs.a0; % Jacobson (2014) radius to equatorial radius
        obs.J2  = +3510.7e-6*cfac^2;
        obs.J4  =   -34.2e-6*cfac^4;
        obs.J6  =          0*cfac^6;
        obs.J8  =          0*cfac^8;
        obs.J10 =          0*cfac^10;
        obs.J12 =          0*cfac^12;
        obs.J14 =          0*cfac^14;

        % Recommended uncertainties
        obs.dJ2  = 0.7e-6*cfac^2;
        obs.dJ4  = 1.3e-6*cfac^4;
        obs.dJ6  = inf;
        obs.dJ8  = inf;
        obs.dJ10 = inf;
        obs.dJ12 = inf;
        obs.dJ14 = inf;

        % It is occasionally convenient to collect the Js and dJs in vectors.
        obs.Js = [-1,obs.J2,obs.J4,obs.J6,obs.J8,...
                     obs.J10,obs.J12,obs.J14];
        obs.dJs = [0,obs.dJ2,obs.dJ4,obs.dJ6,obs.dJ8,...
                     obs.dJ10,obs.dJ12,obs.dJ14];
    end

    function obs = Uranus_ppwd()
        obs = observables.Uranus();
        obs.J6  =  5.1769e-7;
        obs.J8  = -1.0421e-8;
        obs.J10 =  2.5672e-10;
        obs.J12 = -7.2879e-12;
        obs.J14 =  2.4279e-13;
        obs.dJ2  = 1e-6*obs.J2;
        obs.dJ4  = 1e-5*obs.J4;
        obs.dJ6  = 1e-4*obs.J6;
        obs.dJ8  = 1e-4*obs.J8;
        obs.dJ10 = 1e-2*obs.J10;
        obs.dJ12 = 1e-0*obs.J12;
        obs.dJ14 = 1e-0*obs.J14;
        obs.Js = [-1,obs.J2,obs.J4,obs.J6,obs.J8,...
                     obs.J10,obs.J12,obs.J14];
        obs.dJs = [0,obs.dJ2,obs.dJ4,obs.dJ6,obs.dJ8,...
                     obs.dJ10,obs.dJ12,obs.dJ14];
    end

    function obs = Neptune()
        obs = observables();
        obs.pname = 'Neptune';

        % Mass and radius, https://ssd.jpl.nasa.gov/ (2018)
        obs.M = 102.4126e24;
        obs.dM = 0.0048e24;
        obs.a0 = 24764e3;
        obs.s0 = 24622e3;

        % Boundary conditions
        obs.P0 = 1e5;             % The reference radius is the 1 bar level
        obs.T0 = 72; obs.dT0 = 2; % Lindal, G.F., 1992. Astrophys. J. 103
        obs.rho0 = 0.387;         % Ideal gas (mmw=2.319 amu) at (P0,T0)
        obs.drho0 = 0.0108;       % Half the range of T0+/-dT0
        obs.rhomax = 20000;       % Generous guess

        % Rotation rate, https://ssd.jpl.nasa.gov/ (2018)
        obs.P = 0.67125*24*3600;
        obs.w = 2*pi/obs.P;
        obs.GM = obs.G*obs.M;
        obs.q = obs.w^2*obs.a0^3/obs.GM;
        obs.m = obs.w^2*obs.s0^3/obs.GM;
        obs.dP = 600; % Basically a guess, see e.g. Podolak and Helled 2012
        obs.dw = 2*pi/obs.P^2*obs.dP;
        obs.dq = 2*obs.w*obs.a0^3/obs.GM*obs.dw;
        obs.dm = 2*obs.w*obs.s0^3/obs.GM*obs.dw;

        % Some methods (e.g. priors.py) use these fiducial values for scale...
        obs.rhobar = obs.M/(4*pi/3*obs.s0^3);

        % Gravity
        % Nominal coefficients from Jacobson (2009) table 5
        cfac = 25225e3/obs.a0; % Jacobson (2009) radius to equatorial radius
        obs.J2  = +3408.4e-6*cfac^2;
        obs.J4  =   -33.4e-6*cfac^4;
        obs.J6  =          0*cfac^6;
        obs.J8  =          0*cfac^8;
        obs.J10 =          0*cfac^10;
        obs.J12 =          0*cfac^12;
        obs.J14 =          0*cfac^14;

        % Recommended uncertainties
        obs.dJ2  = 4.5e-6*cfac^2;
        obs.dJ4  = 2.9e-6*cfac^4;
        obs.dJ6  = inf;
        obs.dJ8  = inf;
        obs.dJ10 = inf;
        obs.dJ12 = inf;
        obs.dJ14 = inf;

        % It is occasionally convenient to collect the Js and dJs in vectors.
        obs.Js = [-1,obs.J2,obs.J4,obs.J6,obs.J8,...
                     obs.J10,obs.J12,obs.J14];
        obs.dJs = [0,obs.dJ2,obs.dJ4,obs.dJ6,obs.dJ8,...
                     obs.dJ10,obs.dJ12,obs.dJ14];
    end

    function obs = Neptune_ppwd()
        obs = observables.Neptune();
        obs.J6  =  5.7563e-7;
        obs.J8  = -1.2049e-8;
        obs.J10 =  2.9915e-10;
        obs.J12 = -8.2982e-12;
        obs.J14 =  2.6246e-13;
        obs.dJ2  = 1e-6*obs.J2;
        obs.dJ4  = 1e-5*obs.J4;
        obs.dJ6  = 1e-4*obs.J6;
        obs.dJ8  = 1e-4*obs.J8;
        obs.dJ10 = 1e-2*obs.J10;
        obs.dJ12 = 1e-0*obs.J12;
        obs.dJ14 = 1e-0*obs.J14;
        obs.Js = [-1,obs.J2,obs.J4,obs.J6,obs.J8,...
                     obs.J10,obs.J12,obs.J14];
        obs.dJs = [0,obs.dJ2,obs.dJ4,obs.dJ6,obs.dJ8,...
                     obs.dJ10,obs.dJ12,obs.dJ14];
    end
end
end
