################################################################################
# Some functions for looking at sampels of cooked planets
################################################################################

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

def hist_moi(fname, newfig=True, bins='auto', density=True, **kwargs):
    # Prepare the data
    planets = load_planets(fname)
    ice = np.array([p.NMoI for p in planets])

    # Prepare the canvas
    if newfig:
        plt.figure(figsize=(8,6))

    # Plot the histogram
    plt.hist(ice, bins=bins, density=density, **kwargs)

    # Style, annotate, and show
    plt.xlabel(r'Normalized moment of inertia, $I/Ma_0^2$')
    plt.xlim(min(ice),max(ice))
    plt.yticks([])
    plt.show(block=False)
