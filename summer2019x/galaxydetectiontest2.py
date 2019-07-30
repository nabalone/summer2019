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
THRESHOLD = 0.5
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
    
    

        

for filename in filenames[220:221]:
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
    
    objects = sep.extract(swappedData, THRESHOLD, 
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
    flux, _fluxerr, _flag = sep.sum_ellipse(swappedData, objects['x'], 
                                            objects['y'], objects['a'], 
                                            objects['b'], objects['theta'], 
                                            2.5*kronrad, subpix=1)
    magnitude = -2.5 * np.log10(flux/float(image_file[0].header['MJD-OBS'])) + MAG_ZERO 
    
    #size is half light radius
    size, _flag = sep.flux_radius(swappedData, objects['x'], objects['y'], 
                      6.*objects['a'], 0.5, normflux=flux, subpix=5)
        
    separation = np.sqrt(((objects['x']-120)**2) + ((objects['y']-120)**2))
    
    # Observed number density of galaxies brighter than magnitude M (From Berger 2010)
    sigma = 10 ** (0.33 * (magnitude - 24) - 2.44) / (0.33 * np.log(10))

    # Effective radius
    R_effective = np.sqrt(np.abs(separation) + 4 * np.abs(size) ** 2)

    # Probability of chance coincidence
    chanceCoincidence = 1 - np.exp(-np.pi * R_effective ** 2 * sigma)
#TODO check reasonable chance coincidence 
    for i in range(len(objects)):
        if separation[i] > 25:
            chanceCoincidence[i] = 2
    bestCandidate = np.argmin(chanceCoincidence)
    image_file.close()
    

    '''
        if separation[bestCandidate] > 20:
            #print("far: " + filename + " "+ str(separation[bestCandidate]))
            pass
        ellipticity = 1- (objects['b'][bestCandidate]/objects['a'][bestCandidate])
        
    
    
    #TODO make sure directions are correct
        centerRa = image_file[0].header['CRVAL1']
        centerDec = image_file[0].header['CRVAL2']
        ra = centerRa-image_file[0].header['CDELT1'] * objects['x'][bestCandidate]
        dec = centerDec-image_file[0].header['CDELT2']*objects['y'][bestCandidate]
        return [kronrad[bestCandidate], objects['x'][bestCandidate], objects['y'][bestCandidate],
                ra, dec, magnitude[bestCandidate],
                objects['theta'][bestCandidate], ellipticity]
    
    
    #TODO randomize
    
        objects, biggestObjNum = getData(filename)
    '''
    # plot background-subtracted image
    fig, ax = plt.subplots()

    im = ax.imshow(swappedData, interpolation='nearest', cmap='gray')
    
    # plot an ellipse for each object
    #i = biggestObjNum
    i=bestCandidate
    e = Ellipse(xy=(objects['x'][i], objects['y'][i]),
                width=6*objects['a'][i],
                height=6*objects['b'][i],
                angle=objects['theta'][i] * 180. / np.pi)
    e.set_facecolor('none')
    e.set_edgecolor('red')
    ax.add_artist(e)
    plt.savefig(namegen(), dpi=150)
    

