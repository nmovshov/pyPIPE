classdef ppwd
%% Methods for working with the PPWD density parameterization
methods(Static)
    function [zvec,dvec] = ppwd_profile(N,x,rho0,zvec,forcemono)
        %PPWD_PROFILE Polynomial plus two smoothed step functions.
        %    [zvec,dvec] = PPWD_PROFILE(N, x, rho0) returns an N-point density
        %     profile where density d is approximated by a single polynomial of
        %     normalized radius z=s/s0 overlain with two sigmoid functions. The
        %     default radii spacing is one of equal increments between z=1 and
        %     z=1/N. The parameter vector x and the surface density rho0
        %     together fully define the density profile.
        %
        %     The first three elements of x define the location, scale, and
        %     sharpness of the first (inner) step. Specify location 0<z_1<1 in
        %     normalized mean radius and specify scale in kg/m^3. The sharpness
        %     parameter is non-dimensional positive; a larger value results in
        %     a narrower step. Try ~100. The next three elements of x specify
        %     the same parameters for the second (outer) step. The remaining
        %     elements are the coefficients of the polynomial. Remember that a
        %     deg-n polynomial will be constrained by boundary conditions and
        %     defined with n-1 free coefficients, rather than n+1.
        %
        %     You can BYO radii distribution with PPWD_PROFILE(..., zvec), in
        %     which case zvec had better be a 1d vector and N will be ignored.
        %     The returned profile is defined on normalized radii, even if zvec
        %     wasn't.
        %
        %     PPWD_PROFILE(..., forcemono=True) forces the resulting density
        %     profile to be monotonically non increasing. This shouldn't be
        %     necessary but occasionally for high-resolution models and
        %     high-order polynomials there is an unwanted local min somewhere
        %     in the first few layers.
        arguments
            N (1,1) double
            x (1,:) double
            rho0 (1,1) double {mustBeNonnegative}
            zvec (:,1) double {mustBePositive} = []
            forcemono (1,1) logical = true
        end
        y = ppwd.fixpoly(x(7:end), rho0);
        dprof = ppwd.polynomial_profile(N, y, zvec, forcemono);
        ro0 = dprof(1,2);
        dprof = ppwd.add_density_jump(dprof, x(4), x(5), x(6)); % outer
        dprof = ppwd.add_density_jump(dprof, x(1), x(2), x(3)); % inner
        dvec = dprof(:,2) - dprof(1,2) + ro0;
        zvec = dprof(:,1);
        % Undocumented convenience hack: return as 1 array if nargout<2
        if nargout < 2, zvec = [zvec, dvec]; end
    end

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

    function lp = ppwd_prior(x,obs)
        % A prior on the PPWD sampling-space parameters.
        %
        % We use uniform priors on the physical values and we transform those
        % to the correct prior on the sampling-space values. If
        %     Y = log(X - a) - log(b - X)
        % then
        %     X = (a + b*e^Y)/(1 + e^Y)
        % and
        %     X ~ U(a,b)
        % becomes
        %     Y ~ (b - a)*(e^y)/(1 + e^y)^2

        if nargin < 2, obs = []; end
        lp = zeros(size(x));
        s = ppwd.supports(obs);

        % inner z
        a = s.z1(1); b = s.z1(2);
        lp(1) = x(1) - 2*log(1 + exp(x(1))) + log(b - a);

        % inner drho
        a = s.dro1(1); b = s.dro1(2);
        lp(2) = x(2) - 2*log(1 + exp(x(2))) + log(b - a);

        % inner s
        a = s.s1(1); b = s.s1(2);
        lp(3) = x(3) - 2*log(1 + exp(x(3))) + log(b - a);

        % outer z
        a = s.z2(1); b = s.z2(2);
        lp(4) = x(4) - 2*log(1 + exp(x(4))) + log(b - a);

        % outer drho
        a = s.dro2(1); b = s.dro2(2);
        lp(5) = x(5) - 2*log(1 + exp(x(5))) + log(b - a);

        % outer s
        a = s.s2(1); b = s.s2(2);
        lp(6) = x(6) - 2*log(1 + exp(x(6))) + log(b - a);

        % polynomial coefficients
        a = s.poly(1); b = s.poly(2);
        for k=7:length(x)
            lp(k) = x(k) - 2*log(1 + exp(x(k))) + log(b - a);
        end

        % Sum independent priors
        lp = sum(lp);
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
        %   with different values to see the effect. Default 200.
        arguments
            dprof (:,2) double
            z (1,1) double
            scale (1,1) double
            sharpness (1,1) double = 200
        end
        x = dprof(:,1)/dprof(1); % make sure it's normalized
        y = ppwd.sigmoidt(x,z,scale,sharpness);
        dprof(:,2) = dprof(:,2) + y;
    end

    function y = logit(x, sup)
        a = sup(1); b = sup(2);
        y = log(x - a) - log(b - x);
    end

    function y = expit(x, sup)
        a = sup(1); b = sup(2);
        y = (a + b.*exp(x))./(1 + exp(x));
    end

    function y = sigmoidt(x, z, scale, sh)
        y = scale*(pi/2 + atan(-sh*(x - z)))/pi;
    end

    function y = sigmoide(x, z, scale, sh)
        y = scale./(1 + exp(sh*(x - z)));
    end

    function y = ppwd_transform(x, obs)
        if nargin < 2, obs = []; end
        s = ppwd.supports(obs);
        y1 = [ppwd.expit(x(1),s.z1),...
              ppwd.expit(x(2),s.dro1),...
              ppwd.expit(x(3),s.s1)];
        y2 = [ppwd.expit(x(4),s.z2),...
              ppwd.expit(x(5),s.dro2),...
              ppwd.expit(x(6),s.s2)];
        y3 = ppwd.expit(x(7:end),s.poly);
        y = [y1, y2, y3];
    end

    function y = ppwd_untransform(x, obs)
        % Transform ppwd params vector to sample space.
        s = ppwd.supports(obs);
        y1 = [ppwd.logit(x(1),s.z1),...
              ppwd.logit(x(2),s.dro1),...
              ppwd.logit(x(3),s.s1)];
        y2 = [ppwd.logit(x(4),s.z2),...
              ppwd.logit(x(5),s.dro2),...
              ppwd.logit(x(6),s.s2)];
        y3 = ppwd.logit(x(7:end),s.poly);
        y = [y1, y2, y3];
    end

    function a = quadratic_planet(M, R, rho0)
        % Return a s.t. rho(r) = a*(r/R)^2 - a integrates to M.
        a = 5/2*rho0 - 15*M/8/pi/R^3;
    end

    function s = supports(obs)
        if nargin < 1, obs = []; end
        s.z1 = [0.05, 0.5];
        s.z2 = [0., 0.85];
        s.dro1 = [0, 3e4];
        s.dro2 = [0, 3e4];
        s.s1 = [20,1001];
        s.s2 = [20,1001];
        s.poly = [-1e7, 1e7];
        if ~isempty(obs)
            s.dro1 = [0,obs.rhomax];
            s.dro2 = [0,obs.rhomax];
        end
    end
end % methods
end % classdef
