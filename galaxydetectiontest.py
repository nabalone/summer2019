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
#TODO add to csvwriter
import numpy as np
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)

from astropy.io import fits

import glob, sys, os
import csv
import sep
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)

DESTDIR = os.getcwd() 
SOURCEDIR = "F:\ps1hosts"
#? sys.path.append(os.path.abspath(SOURCEDIR))
START_THRESH = 1 # threshhold brightness for detection
THRESH_STEP = 0.3 #amount to decrease threshhold if nothing found
THRESH_MIN = 0.1 #threshhold at which to give up and return filler
MINAREA = 3 #minimum object area for detection
FILLER = [None]*8
AREA_COMP_FACTOR = 1.5
PEAK_COMP_FACTOR = 1.5
#TODO add comparison?

filenames = glob.glob(SOURCEDIR + '\\psc*.fits')


namecount = 0
def namegen():
    global namecount
    namecount += 1
    return DESTDIR + "\galaxyimage" + str(namecount) + ".png"
    
    

#TODO randomize
for filename in filenames[204:207]:
        
    image_file = fits.open(filename)
    image_data = image_file[0].data
    
    #fix byte order
    swappedData = image_data.byteswap(inplace=True).newbyteorder()

  
    # subtracting out background
    bkg = sep.Background(swappedData)
    swappedData = swappedData - bkg 
    #also add as param to extract call below
    
    
    
    #helper function to extract the data from a file
    def getData(filename, threshold=START_THRESH):
        # no object detected
        if threshold < 0.1:
            return FILLER
        objects = sep.extract(swappedData, threshold, err=bkg.globalrms, minarea = MINAREA)
 
#TODO unnecessary?
        # if no sources detected, decrease threshold and try again. If threshold 
        # reaches threshhold min, send filler data
        if len(objects['x']) == 0: #there are no objects
            curthresh = threshold
            return getData(filename, threshhold=curthresh-TRESH_STEP)
        
        #TODO necessary? if so move import statement
        biggestObjNum = np.argmax(objects['npix']) #object with the most pixels
        print(objects['x'][biggestObjNum])
        print(objects['y'][biggestObjNum])
        print('')
        return (objects, biggestObjNum)

        # set low threshhold, low minsize, pick biggest/brightest/centermost with 
        # some weighting. should we ignore brightness? what if its tiny?
        
    objects, biggestObjNum = getData(filename)
    
    # plot background-subtracted image
    fig, ax = plt.subplots()

    im = ax.imshow(swappedData, interpolation='nearest', cmap='gray')
    
    # plot an ellipse for each object
    #i = biggestObjNum
    i=biggestObjNum
    e = Ellipse(xy=(objects['x'][i], objects['y'][i]),
                width=6*objects['a'][i],
                height=6*objects['b'][i],
                angle=objects['theta'][i] * 180. / np.pi)
    e.set_facecolor('none')
    e.set_edgecolor('red')
    ax.add_artist(e)
    plt.savefig(namegen(), dpi=150)

