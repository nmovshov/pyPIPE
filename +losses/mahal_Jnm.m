function d = mahal_Jnm(mdl, obs, degs)
%MAHAL_JNM Weighted Mahalanobis distance in J space.
%   d = MAHAL_JNM(mdl, obs, ords) returns the Mahalanobis distance between the
%   model gravity in mdl and the observed gravity in obs. The mdl parameter is
%   probably an instance of TOFPlanet or CMSPlanet, but all that is required is
%   that is has a property/field Jn for every degree n specified in the vector
%   degs (default degs=[2,4,6]) The paramater obs is probably a struct output from
%   the +observables package, but all that is required is the same fields as mdl
%   and in addition a cov field holding a covariance matrix. The n-by-n obs.cov
%   matrix is assumed to hold covariances sigJ_nm for all even degrees up to 2n.
%   The input vector degs specifies the degrees to include in the distance
%   calculation, the default is degrees 2, 4, and 6, or degs=[2,4,6]. For example,
%   if degs=[2,4] then
%
%       x = [mdl.J2, mdl.J4];
%       mu = [obs.J2, obs.J4];
%   and
%       d = sqrt([x - mu]*inv(mdl.cov)*[x - mu]');
%
%   This loss function is suitable for use when the J uncertainties are
%   correlated. A log-likelihood function proportional to a multivariate normal
%   can simply return negative one half times the square of this loss function.

if nargin < 2
    fprintf('Usage:\n\tmahal_Jnm(mdl, obs, <degs=[2,4,6]>)\n')
    return
end
if ~isfield(obs,'cov')
    error('Missing required field obs.cov; if J errors are uncorrelated use losses.euclid_Jnm instead.')
end
if nargin < 3
    degs = [2,4,6];
end
narginchk(2,3)

x = [];
mu = [];
for k=1:length(degs)
    J = ['J', num2str(degs(k))];
    x = [x, mdl.(J)];
    mu = [mu, obs.(J)];
end
inds = degs/2; % We assume the original covariance is from 2 to n
cov = obs.cov(inds,inds);
d = sqrt((x - mu)*inv(cov)*(x - mu)');

end
