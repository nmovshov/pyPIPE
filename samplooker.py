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
