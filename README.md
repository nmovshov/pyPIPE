# pyPIPE
Drivers and analysis helpers for planetary interiors posteriors exploration

This repository holds ready-to-use drivers for MCMC sampling of interior
density profiles, including generic priors and loss functions, observables, and
example seeds and launch scripts. Typical usage is to clone or export the repo
to a project directory and work from there. A particular project will require
making customized seeds, choosing a sampling strategy, and interpreting the
samples. Very likely each project will tweak one or more of the needed
components: the model itself, or the prior, or the driver, or maybe the loss
function or something else. For help with analysis the ahelpers module
implements many generic tasks (reading a saved chain, plotting traces and
other diagnostics, and more) and will very likely be tweaked as well.

The repo contains local copies of cms and tof solvers from CMS-planet and
TOF-planet. These are not submodules; they are periodically synced with the
latest versions. When in doubt you can always get a fresh copy.

There is an example slurm launch script that can be used on lux.
