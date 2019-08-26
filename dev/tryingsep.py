# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 08:44:34 2019

@author: Faith
"""

import numpy as np
import sep

from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['figure.figsize'] = [10., 8.]
data = fits.open("psc470047.3.fits")
data0 = data[0]
m, s = np.mean(data), np.std(data)
plt.imshow(data, interpolation='nearest', cmap='gray', vmin=m-s, vmax=m+s, origin='lower')
plt.colorbar();