The samples in this folder use the ppwd density profile parameterization with
*rotation non-preserving* gravity solution. The rotation parameter m is fixed
and the implied rotation period is unpredictable, depending on what the
equilibrium shape turns out to be. In all samples
the equatorial radius and 1-bar density are fixed at their reference values and
the implied mass is part of the likelihood. The chi2 dof are then the
number of J values considered plus 1. The sample file naming scheme is

`<poly_deg>_<toforder>_<J_max>_<obsnum>.txt`

This is an attempt to generate models similar to those found with the
hierarchical sampling of piecewise-quadratic and quartic top we used in 2019.
But I don't think this will work, because it was more likely the hierarchical
fixed-jumps nature and/or the different tempering algorithm that were
responsible.