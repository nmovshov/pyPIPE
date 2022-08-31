###############################################################################
# Some functions for looking at sampels of cooked planets
###############################################################################

import sys, os
import pickle
import numpy as np
import matplotlib.pyplot as plt
import ppwd
import TOFPlanet

def load_planets(fname):
    with open(fname, 'rb') as f:
        planets = pickle.load(f)
        print(f"Found {len(planets)} planets in {fname}.")
    return planets

def hist_mass(fname, newfig=True, bins='auto', density=True, **kwargs):
    # Prepare the data
    planets = load_planets(fname)
    M = np.array([p.M for p in planets])/1e24

    # Prepare the canvas
    if newfig:
        plt.figure(figsize=(8,6))
    else:
        plt.figure(plt.gcf().number)

    # Plot the histogram
    plt.hist(M, bins=bins, density=density, **kwargs)

    # Style, annotate, and show
    plt.xlabel(r'Integrated mass $M$ [$10^{24}$ kg]')
    # plt.xlim(M.mean()-2*M.std(),M.mean()+2*M.std())
    plt.yticks([])
    plt.show(block=False)

def hist_moi(fname, newfig=True, bins='auto', density=True, **kwargs):
    # Prepare the data
    planets = load_planets(fname)
    ice = np.array([p.NMoI for p in planets])

    # Prepare the canvas
    if newfig:
        plt.figure(figsize=(8,6))
    else:
        plt.figure(plt.gcf().number)

    # Plot the histogram
    plt.hist(ice, bins=bins, density=density, **kwargs)

    # Style, annotate, and show
    plt.xlabel(r'Normalized moment of inertia, $I/Ma_0^2$')
    # plt.xlim(ice.mean()-2*ice.std(),ice.mean()+2*ice.std())
    plt.yticks([])
    plt.show(block=False)

def hist_J2(fname, newfig=True, bins='auto', density=True, **kwargs):
    # Prepare the data
    planets = load_planets(fname)
    J2 = 1e6*np.array([p.Js[1] for p in planets])

    # Prepare the canvas
    if newfig:
        plt.figure(figsize=(8,6))

    # Plot the histogram
    plt.hist(J2, bins=bins, density=density, **kwargs)

    # Style, annotate, and show
    plt.xlabel(r'$J_2\times{10^6}$')
    # plt.xlim(J2.mean()-2*J2.std(), J2.mean()+2*J2.std())
    plt.yticks([])
    plt.show(block=False)

def ensemble_of_profs(fname, newfig=True, nlines=20, alfa=0.4, **kwargs):
    # Prepare the data
    planets = load_planets(fname)
    profs = np.array([p.rhoi for p in planets]).T
    rcs = profs[-1,:]
    ind = np.argsort(rcs)
    profs = profs[:,ind]

    # Prepare the canvas
    if newfig:
        plt.figure(figsize=(8,6))

    # Plot the lines
    x = planets[0].si/planets[0].s0
    x = np.append(x, 0)
    skip = int(np.ceil(len(planets)/nlines))
    for k in range(0,len(planets),skip):
        y = profs[:,k]
        y = np.append(y,y[-1])
        plt.plot(x,y,alpha=alfa,**kwargs)

    # Style, annotate, and show
    plt.xlabel(r'Level surface radius, $s/R_m$')
    plt.ylabel(r'$\rho$ [1000 kg/m$^3$]')
    plt.show(block=False)

def density_envelope(fname, newfig=True, prctile=2, **kwargs):
    # Prepare the data
    planets = load_planets(fname)
    profs = np.array([p.rhoi for p in planets]).T
    prcs_lo = prctile
    prcs_hi = 100 - prcs_lo
    x = planets[0].si/planets[0].s0
    ylo = np.percentile(profs, prcs_lo, axis=1)/1000
    yhi = np.percentile(profs, prcs_hi, axis=1)/1000

    # Prepare the canvas
    if newfig:
        plt.figure(figsize=(8,6))

    # Plot the shaded regions
    plt.fill_between(x, ylo, yhi, **kwargs)

    # Style, annotate, and show
    plt.xlabel(r'Level surface radius, $s/R_m$')
    plt.ylabel(r'$\rho$ [1000 kg/m$^3$]')
    plt.show(block=False)

def barotrope_envelope(fname, newfig=True, prctile=2, **kwargs):
    from scipy.interpolate import interp1d
    # Prepare the data
    planets = load_planets(fname)
    pees = 1e-11*np.array([p.Pi for p in planets]).T
    rhos = 1e-3*np.array([p.rhoi for p in planets]).T
    x = np.logspace(-6, np.ceil(np.log10(pees.max())), 1024)
    y = np.nan*np.ones((1024, len(planets)))
    for k in range(len(planets)):
        rhoofp = interp1d(pees[:,k], rhos[:,k], kind='cubic')
        ind = x < pees[-1,k]
        y[ind,k] = rhoofp(x[ind])

    prcs_lo = prctile
    prcs_hi = 100 - prcs_lo
    ylo = np.nanpercentile(y, prcs_lo, axis=1)
    yhi = np.nanpercentile(y, prcs_hi, axis=1)

    # Prepare the canvas
    if newfig:
        plt.figure(figsize=(8,6))

    # Plot the shaded regions
    ind = ~np.isnan(ylo)
    plt.fill_between(x[ind], ylo[ind], yhi[ind], **kwargs)

    # Style, annotate, and show
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(left=1e-5)
    plt.ylim(1e-3, 100)
    plt.vlines(pees[-1,:].min(), *plt.ylim(), linestyle=':', color='r')
    plt.xlabel(r'$p$ [Mbar]')
    plt.ylabel(r'$\rho$ [1000 kg/m$^3$]')
    plt.show(block=False)

def plot_profile(s, rho, newfig=True, **kwargs):
    # Prepare the data
    x = np.append(s, 0)/s[0]
    y = np.append(rho, rho[-1])/1000

    # Prepare the canvas
    if newfig:
        plt.figure(figsize=(8,6))

    # Plot the lines
    plt.plot(x, y, **kwargs)

    # Style, annotate, and show
    plt.xlabel(r'Level surface radius, $s/R_m$')
    plt.ylabel(r'$\rho$ [1000 kg/m$^3$]')
    plt.show(block=False)

def hist_mass_err(fname,obs,newfig=True,bins='auto',density=True,**kwargs):
    # Prepare the data
    planets = load_planets(fname)
    M_err = (np.array([p.M for p in planets]) - obs.M)/obs.dM

    # Prepare the canvas
    if newfig:
        plt.figure(figsize=(8,6))
    else:
        plt.figure(plt.gcf().number)

    # Plot the histogram
    plt.hist(M_err, bins=bins, density=density, **kwargs)

    # Style, annotate, and show
    plt.xlabel(r'$(M-M^\mathrm{obs})/\sigma{M}^\mathrm{obs}$')
    # plt.xlim(M.mean()-2*M.std(),M.mean()+2*M.std())
    plt.yticks([])
    plt.show(block=False)

def hist_J2_err(fname,obs,newfig=True,bins='auto',density=True,**kwargs):
    # Prepare the data
    planets = load_planets(fname)
    J_err = (np.array([p.Js[1] for p in planets]) - obs.J2)/obs.dJ2

    # Prepare the canvas
    if newfig:
        plt.figure(figsize=(8,6))
    else:
        plt.figure(plt.gcf().number)

    # Plot the histogram
    plt.hist(J_err, bins=bins, density=density, **kwargs)

    # Style, annotate, and show
    plt.xlabel(r'$(J_2-J_2^\mathrm{obs})/\sigma{J_2}^\mathrm{obs}$')
    # plt.xlim(M.mean()-2*M.std(),M.mean()+2*M.std())
    plt.yticks([])
    plt.show(block=False)
