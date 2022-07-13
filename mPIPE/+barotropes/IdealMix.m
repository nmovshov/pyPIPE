classdef IdealMix < barotropes.Barotrope
    %IDEALMIX An ideal mixture of N other barotropes.
    
    %% Properties
    properties (Access = public)
        EOSs % vector of barotrope.Barotrope objects to mix
        mrats % vector of mass ratios
    end
        
    %% The constructor
    methods
        function obj = IdealMix(eoss, xs)
            if nargin > 0
                validateattributes(eoss,{'barotropes.Barotrope'},{'vector'})
                validateattributes(xs,{'numeric'},{'nonnegative','finite','vector'})
                assert(length(eoss) == length(xs))
                assert(abs(sum(xs) - 1) < eps)
                obj.EOSs = eoss;
                obj.mrats = xs;
            end
        end
    end
    
    %% Required barotrope methods
    methods
        function PF = test(~)
            PF = true;
        end
        
        function P = pressure(obj,rho)
            % Mix pressure by Dalton's law.
            assert(isvector(rho));
            P = zeros(size(rho));
            for k=1:length(obj.mrats)
                P = P + obj.EOSs(k).pressure(rho)*obj.mrats(k);
            end
        end
        
        function rho = density(obj,P)
            % Mix density by Amagat's law.
            assert(isvector(P));
            invrho = zeros(size(P));
            for k=1:length(obj.mrats)
                invrho = invrho + obj.mrats(k)./obj.EOSs(k).density(P);
            end
            rho = 1./invrho;
        end
    end
    
end
