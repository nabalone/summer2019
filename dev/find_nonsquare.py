# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 14:27:04 2019

@author: Faith
"""

import glob    
from astropy.io import fits

filenames = sorted(glob.glob(SOURCEDIR + '/psc*.6.fits'))
x = []
for filename in filenames[:200]:
    with fits.open(filename) as image_file:
    
        if image_file[0].header['NAXIS1'] != 240 or image_file[0].header['NAXIS2'] != 240:
            x.append(filename)
print x