# Set workspace (paths and/or common variables).

import sys, os
import numpy as np
if sys.platform == 'win32':
    import matplotlib.pyplot as plt

sys.path.append(os.getcwd())

import observables
import ppwd
import TOFPlanet
import ahelpers as ah
import losses
if sys.platform == 'win32':
    import samplooker as spl
