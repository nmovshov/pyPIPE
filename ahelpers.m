classdef ahelpers
    %AHELPERS Some helper functions.

methods (Static)

    function M = mass_integral(svec, dvec)
        % Return approximate mass integral.
        M = -4*pi*trapz(svec, dvec.*svec.^2);
    end
    
    function L = lossify_planets(planets,obs,jmax,rho0,mass)
        narginchk(3,5)
        if nargin < 4 | isempty(rho0), rho0 = true; end
        if nargin < 5 | isempty(mass), mass = false; end
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
