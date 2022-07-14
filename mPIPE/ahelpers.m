classdef ahelpers
    %AHELPERS Some helper functions.

methods (Static)
    function M = mass_integral(svec, dvec)
        % Return approximate mass integral.
        M = -4*pi*trapz(svec, dvec.*svec.^2);
    end
end
end
