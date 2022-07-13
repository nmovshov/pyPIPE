function y = rhomax(cmp, obs)
% Penalty for unreasonable central density, to steer fminsearch.

romax = cmp.rhoi(end);
y = (max(romax, obs.rhomax) - obs.rhomax)/obs.rhobar;
end
