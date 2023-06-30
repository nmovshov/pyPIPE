classdef ppbs
%% Methods for working with the PPBS barotrope parameterization
methods(Static)
    function p = ppbs_planet(N,x,obs)
        %PPWD_PLANET Create a TOFPlanet object from sample-space parameters.
        %   p = PPWD_PLANET(N,x,obs) returns a TOFPlanet object initialized
        %   with observables obs and with N-level density ppwd_profile derived
        %   from sample-space parameters x.
        y = ppwd.ppwd_transform(x,obs);
        [zvec,dvec] = ppwd.ppwd_profile(N,y,obs.rho0);
        p = TOFPlanet(obs);
        p.si = zvec*obs.s0;
        p.rhoi = dvec;
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

    function y = logit(x, sup)
        a = sup(1); b = sup(2);
        y = log(x - a) - log(b - x);
    end

    function y = expit(x, sup)
        a = sup(1); b = sup(2);
        y = (a + b.*exp(x))./(1 + exp(x));
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
        if nargin < 3, rho0 = 0; end
        a = 5/2*rho0 - 15*M/8/pi/R^3;
    end

    function s = supports(obs)
        if nargin < 1, obs = []; end
        s.z1 = [0.05, 0.5];
        s.z2 = [0.5, 0.85];
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
