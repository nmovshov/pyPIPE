function y = mass(cmp, obs)
% sqrt((M_calc - M_obs)^2/dM_obs^2)

narginchk(2,2)

y = sqrt(((cmp.M - obs.M)/obs.dM)^2);
end
