# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 06:55:23 2019

@author: Faith
"""

# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-

# http://docs.astropy.org/en/stable/generated/examples/io/plot_fits-image.html
# http://docs.astropy.org/en/stable/io/fits/

"""
=======================================
Read and plot an image from a FITS file
=======================================

This example opens an image stored in a FITS file and displays it to the screen.

This example uses `astropy.utils.data` to download the file, `astropy.io.fits` to open
the file, and `matplotlib.pyplot` to display the image.

-------------------

*By: Lia R. Corrales, Adrian Price-Whelan, Kelle Cruz*

*License: BSD*

-------------------


"""
import numpy as np
from matplotlib.patches import Ellipse
import os
import glob
import csv
import sep
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits

SOURCEDIR = "C:/Users/Faith/Desktop/noey2019summer/ps1hosts"
DESTDIR = os.getcwd() + "\galaxiestestpix"
FILLER_VAL = None
THRESHOLD = 1.5
MINAREA = 5
MAG_ZERO =25 # for calculating magnitude from flux


# make header
HEADER =['ID']
perImageHeaders = ['KronRad', 'x', 'y', 'RA', 'DEC', 'KronMag', 'Angle', 'Ellipticity']
for i in range(3,7):
    for val in perImageHeaders:
        HEADER.append(val + '_' + str(i))


filenames = sorted(glob.glob('C:/Users/Faith/Desktop/noey2019summer/ps1hosts/psc*'))

namecount = 0
def namegen():
    global namecount
    namecount += 1
    return DESTDIR + "\galaxyimage" + str(namecount) + ".png"
    
    
messfiles = ['C:/Users/Faith/Desktop/noey2019summer/ps1hosts\\psc000010.3.fits', 
             'C:/Users/Faith/Desktop/noey2019summer/ps1hosts\\psc000010.4.fits', 
             'C:/Users/Faith/Desktop/noey2019summer/ps1hosts\\psc000010.5.fits', 
             'C:/Users/Faith/Desktop/noey2019summer/ps1hosts\\psc000010.6.fits', 
             'C:/Users/Faith/Desktop/noey2019summer/ps1hosts\\psc000190.3.fits', 
             'C:/Users/Faith/Desktop/noey2019summer/ps1hosts\\psc000190.4.fits', 
             'C:/Users/Faith/Desktop/noey2019summer/ps1hosts\\psc000190.5.fits',
             'C:/Users/Faith/Desktop/noey2019summer/ps1hosts\\psc000190.6.fits',
             'C:/Users/Faith/Desktop/noey2019summer/ps1hosts\\psc010163.3.fits',
             'C:/Users/Faith/Desktop/noey2019summer/ps1hosts\\psc010163.4.fits',
             'C:/Users/Faith/Desktop/noey2019summer/ps1hosts\\psc010163.5.fits', 
             'C:/Users/Faith/Desktop/noey2019summer/ps1hosts\\psc010163.6.fits']
        
# making blocking array for 010163
arr = []
for i in range(40):
    row = []
    for j in range(240):
        row.append(True)
    arr.append(row)
for i in range(200):
    row = []
    for j in range(240):
        row.append(False)
    arr.append(row)
maskarr = np.array(arr)

for filename in ['C:/Users/Faith/Desktop/noey2019summer/ps1hosts\psc370330.6.fits']:
   # print filename
#TODO does this need to be recursive? if so move the following out, if not maybe
    #get rid of helper fxn
    image_file = fits.open(filename)
    image_data = image_file[0].data

    # https://www.aanda.org/index.php?option=com_article&access=doi&doi=
    # 10.1051/0004-6361/201015362&Itemid=129#FD14
    
    #fix byte order
    swappedData = image_data.byteswap(inplace=True).newbyteorder()
    
    # subtracting out background
    bkg = sep.Background(swappedData)
    swappedData = swappedData - bkg 
    
    objects = sep.extract(swappedData, THRESHOLD, deblend_cont = 1.0, mask = maskarr,
                          err=bkg.globalrms, minarea = MINAREA)
    
    # if no sources detected, decrease threshold and try again. If threshold 
    # reaches threshhold min, send filler data
 #TODO if not objects:
    
    
    '''Find most likely host galaxy'''
    # how to calculate kron radius and flux from 
    # https://sep.readthedocs.io/en/v1.0.x/apertures.html
    kronrad, _krflag = sep.kron_radius(swappedData, objects['x'], objects['y'], 
                                       objects['a'], objects['b'], 
                                       objects['theta'], 6.0)
    print kronrad
    
    
   
    for i in range(len(objects)):
        ki = kronrad[i]
        #check if ki is nan:
        if (not ki == 0) and (not ki < 0) and (not ki > 0):
            print filename
            print objects['peak'][i]

        # plot background-subtracted image
            fig, ax = plt.subplots()
        
            im = ax.imshow(swappedData, interpolation='nearest', cmap='gray')
            
            e = Ellipse(xy=(objects['x'][i], objects['y'][i]),
                        width=6*objects['a'][i],
                        height=6*objects['b'][i],
                        angle=objects['theta'][i] * 180. / np.pi)
            e.set_facecolor('none')
            e.set_edgecolor('red')
            ax.add_artist(e)
            plt.savefig(namegen(), dpi=150)
        

