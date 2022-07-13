classdef ppwd
%% Methods for working with the PPWD density parameterization
methods(Static)
    function [zvec,dvec] = polynomial_profile(N,x,zvec,forcemono)
        %POLYNOMIAL_PROFILE A single polynomial in normalized radius.
        %    [zvec,dvec] = POLYNOMIAL_PROFILE(N, x) returns an N-point density
        %    profile defined by a single polynomial of normalized radius z,
        %    with coefficients x. The default radii spacing is one of equal
        %    increments between z=1 and z=1/N.
        %
        %    You can BYO radii distribution with POLYNOMIAL_PROFILE(N,x,zvec),
        %    in which case zvec had better be a 1d vector and N will be
        %    ignored. The returned profile is defined on normalized radii, even
        %    if zvec wasn't.
        %
        %     POLYNOMIAL_PROFILE(..., forcemono=True) forces the resulting
        %     density profile to be monotonically non increasing. This is
        %     actually the default. It helps when a high degree polynomial
        %     tries hard to asymptote horizontally as x->1 and ends up with a
        %     tiny upward slope.
        arguments
            N (1,1) double
            x (1,:) double
            zvec (:,1) double {mustBePositive} = []
            forcemono (1,1) logical = true
        end

        if isempty(zvec), zvec = linspace(1,1/N,N); end
        dvec = polyval(x,zvec);
        if forcemono
            dvec(1) = max(dvec(1), 0);
            for k=2:N
                dvec(k) = max(dvec(k), dvec(k-1));
            end
        end
        dvec = dvec(:);
        zvec = zvec(:);
        % Undocumented convenience hack: return as 1 array if nargout<2
        if nargout < 2, zvec = [zvec, dvec]; end
    end

    function y = fixpoly(x, rho0, z0, rho0p)
        %FIXPOLY Return full polynomial coefficients from constrained ones.
        %
        %    y = fixpoly(x,rho0) interprets the vector x as the n-1
        %    coefficients of a degree-n polynomial that is constrained to have
        %    no linear term and pass through point (1,rho0). The returned
        %    vector will hold the n+1 coefficients of the full polynomial.
        %
        %    y = fixpoly(x,rho0,z0) constrains that polynomial to pass through
        %    point (z0,rho0) instead.
        %
        %    y = fixpoly(x,rho0,z0,rho0p) additionally constrains the
        %    derivative to pass through (z0,rho0p). Note that this means a
        %    vector x of length n defines a polynomial of degree n+2.
        arguments
            x (1,:) double
            rho0 (1,1) double
            z0 (1,1) double = 1.0
            rho0p (1,1) double = nan
        end
        if isnan(rho0p)
            y = [x, 0, (rho0 - polyval([x,0,0], z0))];
        else
            y = [x, 0, 0, 0];
            a2 = (rho0p - polyval(polyder(y), z0))/(2*z0);
            a0 = rho0 - polyval([x,a2,0,0], z0);
            y = [x,a2,0,a0];
        end
    end

    function dprof = add_density_jump(dprof, z, scale, sharpness)
        %ADD_DENSITY_JUMP Add a localized density increase to a profile.
        %
        %   dprof = ADD_DENSITY_JUMP(dprof, z, scale) takes an existing density
        %   profile (N-by-2 array) and adds a smoothed-step-function localized
        %   in the neighborhood of normalized radius z. The scale parameter
        %   controls the overall height of the step (specified in real density
        %   units). In other words, scale*kg/m^3 will be added to the central
        %   density while nothing will be added to the surface density, with
        %   the bulk of the increase happening near z.
        %
        %   ADD_DENSITY_JUMP(...,sharpness) allows you to control the step's
        %   "sharpness", the radial distance over which most of the density
        %   will be added. The sharpness parameter is dimensionless; experiment
        %   with different values to see the effect. Default 100.
        %
        %   Algorithm: we use the inverse tangent function to approximate the
        %   step.
        arguments
            dprof (:,2) double
            z (1,1) double
            scale (1,1) double
            sharpness (1,1) double = 100
        end
        x = dprof(:,1)/dprof(1); % make sure it's normalized
        y = scale*(pi/2 + atan(-sharpness*(x - z)))/pi;
        dprof(:,2) = dprof(:,2) + y;
        end

end % methods
end % classdef
