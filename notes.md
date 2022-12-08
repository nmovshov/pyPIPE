* Can get more J4_obs0 by starting from J4_obs0 seeds, but must run high temp
  for a while to gurantee new ones. Works with (10000,4000,4000,10000) steps
* In the J4_obs0 sample there are 3 models with J6 within 1e-6 of obs0.J6; but
  so far no deg4 models are able to match J6 at Saturn0 levels (<0.1e-6)
* Do not underestimate how undetermined obs.s0 actually is! The message is
  clear: do not use mrot as anything but an intermediate variable internal to
  the ToF calculation. It's not an input to anything physical.