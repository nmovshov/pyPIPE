function d = euclid_Jnm(mdl, obs, degs)
%EUCLID_JNM Weighted euclidean distance in J space.
%   d = EUCLID_JNM(mdl, obs, ords) returns the eucledean distance between the
%   model gravity in mdl and the observed gravity in obs. The mdl parameter is
%   probably an instance of TOFPlanet or CMSPlanet, but all that is required is
%   that is has a property/field Jn for every degree n specified in the vector
%   degs (default degs=[2,4,6]) The paramater obs is probably a struct output from
%   the +observables package, but all that is required is the same fields as mdl
%   and in addition the corresponding dJn fields holding the 1-sigma uncertainties
%   used to define the distance units. The input vector degs specifies the degrees
%   to include in the distance calculation, the default is degrees 2, 4, and 6, or
%   degs=[2,4,6]. For example, if degs=[2,4] then
%
%       d = sqrt((mdl.J2 - obs.J2)/obs.dJ2)^2 + (mdl.J4 - obs.J4)/obs.dJ4)^2)
%
%   This loss function is suitable for use when the J uncertainties are assumed to
%   be uncorrelated. The dJ values in obs are interpreted as 1-sigma values for
%   the purpose of weighing the distance. In other words, this is equivalent to
%   the Mahalanobis distance with a diagonal covariance. A log-likelihood function
%   proportional to a multivariate normal can simply return negative one half
%   times the square of this loss function.

if nargin < 2
    fprintf('Usage:\n\teuclid_Jnm(mdl, obs, <degs=[2,4,6]>)\n')
    return
end
if nargin < 3
    degs = [2,4,6];
end
narginchk(2,3)

d = 0;
for k=1:length(degs)
    J = ['J', num2str(degs(k))];
    dJ = ['dJ', num2str(degs(k))];
    d = d + ((mdl.(J) - obs.(J))/obs.(dJ))^2;
end

d = sqrt(d);

end
