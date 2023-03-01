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

    function k2 = lovek2(zvec, dvec)
    %LOVEK2 Tidal Love number k2 from density profile.
    %
    % Algorithm
    % ---------
    % Sterne 1939 as described in Buhler 2016.

    % Input handling
    if nargin == 0
        fprintf('Usage:\n\tlovek2(zvec, dvec)\n')
        return
    end
    narginchk(2,2);
    validateattributes(zvec,{'numeric'},{'real','finite','positive','vector'})
    validateattributes(dvec,{'numeric'},{'real','finite','nonnegative','vector'})
    zvec = zvec(:);
    dvec = dvec(:);
    assert(isequal(size(zvec),size(dvec)))
    if zvec(2) < zvec(1)
        zvec = flipud(zvec);
        dvec = flipud(dvec);
    end
    zvec = zvec/zvec(end);

    % Climb up the radius with Buhler (2016) eq. 2
    m = dvec(1)*zvec(1)^3; % starting mass
    rhom = m/zvec(1)^3; % starting mean density
    eta = 0;
    for k=1:length(zvec) - 1
        s1 = (6 - 6*(dvec(k)/rhom)*(eta + 1) + eta - eta^2)/zvec(k);
        zhalf = zvec(k) + 0.5*(zvec(k+1) - zvec(k));
        dhalf = dvec(k) + 0.5*(dvec(k+1) - dvec(k));
        mhalf = m + dhalf*(zhalf^3 - zvec(k)^3);
        rhalf = mhalf/zhalf^3;
        ehalf = eta + s1*(zhalf - zvec(k));
        s2 = (6 - 6*(dhalf/rhalf)*(ehalf + 1) + ehalf - ehalf^2)/zhalf;
        eta = eta + s2*(zvec(k+1) - zvec(k));
        m = mhalf + dvec(k+1)*(zvec(k+1)^3 - zhalf^3);
        rhom = m/zvec(k+1)^3;
    end

    % Return
    k2 = (3 - eta)/(2 + eta);
    end
end
end
