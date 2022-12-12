classdef ahelpers
    %AHELPERS Some helper functions.

methods (Static)
    function z = metal2z(x,Y)
        % METAL2Z Assuming H2/He/H2O mixture convert solar units O:H to Z mass fraction.

        o2h_sol = 0.85e-3; % Solar atmosphere O:H (Anders and Grevesse, 1989, Table 2)
        o2h = x*o2h_sol;
        murat = 2/18; % H2 to H2O mmw ratio

        O2H = @(Z)(Z./((1 - Y)*(1 - Z))*murat)./(2 + 2*murat*Z./((1 - Y)*(1 - Z)));
        z = zeros(size(x));
        for k=1:length(x)
            fun = @(Z)O2H(Z) - o2h(k);
            z(k) = fzero(fun,[0,1 - eps]);
        end
    end
    
    function [YLO, YHI] = density_prctiles(profs, prcs)
        %DENSITY_PRCTILES High and low density values of symmetric percentile.

        prcs_lo = prcs;
        prcs_hi = 100 - prcs_lo;

        % Get the data in the above percentiles
        rhos = profs';
        YLO = prctile(rhos, prcs_lo);
        YHI = prctile(rhos, prcs_hi);
        YLO(:,end+1) = YLO(:,end);
        YHI(:,end+1) = YHI(:,end);
    end

    function M = mass_integral(svec, dvec)
        % Return approximate mass integral.
        M = -4*pi*trapz(svec, dvec.*svec.^2);
    end
    
    function L = lossify_planets(planets,obs,jmax,rho0,mass)
        narginchk(3,5)
        if nargin < 4 || isempty(rho0), rho0 = true; end
        if nargin < 5 || isempty(mass), mass = false; end
        L = nan(size(planets));
        for k=1:length(planets)
            p = planets(k);
            el = 0;
            for j=2:(jmax/2 + 1)
                el = el + ((p.Js(j) - obs.Js(j))/obs.dJs(j))^2;
            end
            if rho0
                el = el + ((p.rho0 - obs.rho0)/obs.drho0)^2;
            end
            if mass
                el = el + ((p.M - obs.M)/obs.dM)^2;
            end
            L(k) = sqrt(el);
        end
    end
    
end
end
