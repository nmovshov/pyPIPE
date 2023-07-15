classdef ppbs
%% Methods for working with the PPBS barotrope parameterization
methods(Static)
    function p = ppbs_planet(N, y, obs, toforder, xlevels)
        %PPBS_PLANET Create a TOFPlanet from ppbs parameters.

        if (nargin == 0) && (nargout == 0)
            fprintf('Usage:\n\tppbs_planet(N, y, obs, toforder={4},')
            fprintf(' xlevels={128}\n')
            return
        end
        if nargin < 4 || isempty(toforder), toforder = 4; end
        if nargin < 5 || isempty(xlevels), xlevels = 128; end

        % The ppbs parameterization
        K1 = y(1); n1 = y(2);
        K2 = y(3); n2 = y(4);
        K3 = y(5); n3 = y(6);
        r12 = y(7); r23 = y(8);

        % Make radius grid; snap to layer boundaries
        zvec = linspace(1, 1/N, N);
        [~, tind] = min(abs(zvec - r12));
        [~, cind] = min(abs(zvec - r23));
        zvec(tind) = r12;
        zvec(cind) = r23;

        % Create the eoss
        eos1 = barotropes.Polytrope(K1, n1);
        eos2 = barotropes.Polytrope(K2, n2);
        eos3 = barotropes.Polytrope(K3, n3);

        % Construct TOFPlanet and assign options and eoss
        p = TOFPlanet(obs);
        p.opts.toforder = toforder;
        p.opts.xlevels = xlevels;
        p.eos = [repmat(eos1, tind - 1, 1);...
                 repmat(eos2, cind - tind, 1);...
                 repmat(eos3, N - cind + 1, 1)];

        % Initialize with quadratic density profile
        a = ppbs.quadratic_planet(obs.M, obs.a0);
        p.si = zvec*obs.a0;
        p.rhoi = a*zvec.^2 - a;
    end


    %% Helper functions
    function y = transform(x, sup)
        % Transform mcmc sample space to ppbs params.
        if nargin < 2 || isempty(sup), sup = ppbs.def_sup; end
        y = nan(size(x));
        for k=1:length(x)
            y(k) = ppbs.expit(x(k), sup{k});
        end
    end

    function x = untransform(y, sup)
        % Transform ppbs params vector to sample space.
        if nargin < 2 || isempty(sup), sup = ppbs.def_sup; end
        x = nan(size(y));
        for k=1:length(y)
            x(k) = ppbs.logit(y(k), sup{k});
        end
    end

    function y = logit(x, sup)
        a = sup(1); b = sup(2);
        y = log(x - a) - log(b - x);
    end

    function y = expit(x, sup)
        a = sup(1); b = sup(2);
        y = (a + b.*exp(x))./(1 + exp(x));
    end

    function a = quadratic_planet(M, R, rho0)
        % Return a s.t. rho(r) = a*(r/R)^2 - a integrates to M.
        if nargin < 3, rho0 = 0; end
        a = 5/2*rho0 - 15*M/8/pi/R^3;
    end

end % methods

properties(Constant=true,Access=public)
    def_sup = {[1e4, 8e5],...    % K1
               [0.02, 2.0],...   % n1
               [1e4, 8e5],...    % K2
               [0.02, 2.0],...   % n2
               [1e4, 8e5],...    % K3
               [0.02, 2.0],...   % n3
               [0.45, 0.999],... % r12
               [0.001, 0.8]};    % r23

    mono_y_seed = [589361.3,...
                   1.2993193778,...
                   589361.3,...
                   1.2993193778,...
                   589361.3,...
                   1.2993193778,...
                   0.6,...
                   0.3];

end % properties

end % class
