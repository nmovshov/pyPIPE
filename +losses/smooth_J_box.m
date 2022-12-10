function d = smooth_J_box(mdl, obs, degs)
%SMOOTH_J_BOX Indifference inside a 1-sigma box, euclidean distance outside.
%   d = SMOOTH_J_BOX(mdl, obs, ords) returns the eucledean distance between the
%   model gravity in mdl.Jn and the nearest corner of the 1-sigma box defined by
%   obs.Jn+/-obs.dJn, for n in degs. The mdl parameter is probably an
%   instance of TOFPlanet or CMSPlanet, but all that is required is that is has a
%   property/field Jn for every degree n specified in the vector degs (default
%   degs=[2,4,6]) The paramater obs is probably a struct output from the
%   +observables package, but all that is required is the same fields as mdl and
%   in addition the corresponding dJn fields holding the 1-sigma uncertainties
%   used to define the distance units. The input vector degs specifies the degrees
%   to include in the distance calculation, the default is degrees 2, 4, and 6, or
%   degs=[2,4,6].
%
%   This loss function is useful when the J uncertainties are derived from a
%   separate process as the nomial J values and we don't wish to privilege the
%   central values over anything else in the error range.

if nargin < 2
    fprintf('Usage:\n\tsmooth_J_box(mdl, obs, <degs=[2,4,6]>)\n')
    return
end
if nargin < 3
    degs = [2,4,6];
end
narginchk(2,3)

WD = [];
for k=1:length(degs)
    J = ['J', num2str(degs(k))];
    dJ = ['dJ', num2str(degs(k))];
    if mdl.(J) > (obs.(J) + obs.(dJ))
        WD = [WD, (mdl.(J) - (obs.(J) + obs.(dJ)))/obs.(dJ)];
    elseif mdl.(J) < (obs.(J) - obs.(dJ))
        WD = [WD, (mdl.(J) - (obs.(J) - obs.(dJ)))/obs.(dJ)];
    else
        WD = [WD, 0];
    end
end

d = sqrt(sum(WD.^2));

end
