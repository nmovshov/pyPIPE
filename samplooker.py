################################################################################
# Some functions for looking at sampels of cooked planets
################################################################################

import sys, os
import pickle
import numpy as np
import matplotlib.pyplot as plt

def load_planets(fname):
    with open(fname, 'rb') as f:
        planets = pickle.load(f)
        print(f"Found {len(planets)} planets in {fname}.")
    return planets